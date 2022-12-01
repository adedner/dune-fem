# SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception

###############################################################################
# This paraview reader adds support for 'dune binary format' # files (dbf).
# The file is assumed to be written using 'dune.common.pickle.dump'. It
# therefore consists of two parts (the required jit module source code and
# a pickled list of objects). This list is searched for objects containing
# a 'gridView' attribute - these are all assumed to be grid functions
# over the same grid view and with a 'pointData' attribute.
# If no entry in the list with a 'gridView' attribute is found the first
# entry is assumed to be a grid view and only the grid is plotted.
###############################################################################

import numpy as np
import os,sys,vtk,importlib,glob,json
from importlib.util import spec_from_loader, module_from_spec
from importlib.machinery import SourceFileLoader
from paraview.util.vtkAlgorithm import VTKPythonAlgorithmBase
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

# patch stdout/stderr to include 'isatty' method to make ufl happy
# (vtk provides its own output streams without that method
# ufl fails if 'AttributeError').
# Note that the pvpython streams have no slots so 'setattr' does not work
# and a more complex hack is required:
def patchStream(stream):
    class PatchedStream(stream.__class__):
        def __init__(self):
            self.__impl__ = stream
        def isatty(self):
            return False
        def __getattr__(self, item):
            def tocontainer(func):
                @wraps(func)
                def wrapper(*args, **kwargs):
                    return func(*args, **kwargs)
                return wrapper
            result = getattr(self.__impl__, item)
            if not isinstance(result, PatchedStdStream) and callable(result):
                result = tocontainer(result)
            return result
    return PatchedStream()
sys.stdout = patchStream(sys.stdout)
sys.stderr = patchStream(sys.stderr)

# In older paraview versions there is no way to set the
# virtual environment to use - use a environment variable
# to set it before starting paraview.
# This should be improved to take other usages into account.

# This finds all 'egg-link' files in a given folder structure.
# These correspond to packages installed 'editable' and need to
# be added by hand to the Python search path:
def find_egglinks(directory_name):
    dune_found = []
    for path, subdirs, files in os.walk(directory_name):
        if not path.endswith("site-packages"):
            continue
        dune_found.append(path)
        for name in files:
            if not "dune" in name:
                continue
            ext = os.path.splitext(name)[1]
            if ext == ".egg-link":
                file_path = os.path.join(path,name)
                with open(file_path,"r") as f:
                    dune_found.append(f.read().split()[0])
    return dune_found
# We find an active virtual env by checking if the environment variable
# 'VIRTUAL_ENV' is set - this work at least with activated environments
# setup with 'venv' on Linux:
def setDuneModulePaths():
    try:
        envdir = os.path.realpath(os.environ['VIRTUAL_ENV'])
        dunePaths = find_egglinks(os.path.join(envdir,"lib"))
        sys.path += dunePaths
        if not "DUNE_PY_DIR" in os.environ:
            os.environ["DUNE_PY_DIR"] = os.path.join(envdir,".cache")
        sys.path += os.path.join(os.environ["DUNE_PY_DIR"],"python","dune","generated")
        # print(os.environ["DUNE_PY_DIR"], dunePaths)
    except KeyError:
        # print("no virtual env path found!")
        pass

############################################



# Actual reader
# -------------
# Some documentation
# https://kitware.github.io/paraview-docs/latest/python/paraview.util.vtkAlgorithm.html
# https://github.com/Kitware/ParaView/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py

dune_extensions = ["dbf","tsdbf"] # ,"dgf"]
@smproxy.reader(
    label="Dune Reader",
    extensions=dune_extensions,
    file_description="dune binary format files",
)
class DuneReader(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid"
        )
        self._filename = None
        self._filenameSeries = None
        self._level = 0
        self._transform = None
        self._transformFcts = []
        self._transformFct = ""
        self._dataFcts = []
        self._dataFct = 0
        self._timeSteps = None
        self._currentTime = 0
        self._gridView = None
        setDuneModulePaths()
        try:
            import dune.common.pickle
            import dune.common.utility
        except ImportError:
            raise ImportError("could not import dune.common")
        self.load = dune.common.pickle.load
        self.reload = dune.common.utility.reload_module

    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser( extensions=dune_extensions, file_description="dune binary file format" )
    def SetFileName(self, filename):
        if (self._filename != filename):
            self._filename = filename
            if self._filename != "None":
                filepart = filename.split(".")
                if len(filepart)>=3 and filepart[2] == "dbf":
                    if filepart[1] == "series":
                        with open(filename,"r") as f:
                            self._filenameSeries = json.load(f)
                        print(self._filenameSeries)
                        self._timeSteps = [float(v["time"]) for v in self._filenameSeries.values()]
                        self._currentTime = self._timeSteps[0]
                    else:
                        # assume a file of the form 'base.0000.dbf' was
                        # provided and find all files in the numbered series
                        self._filenameSeries = glob.glob(".".join(filepart[0:-2]) + ".*.dbf")
                        self._filenameSeries.sort()
                        self._timeSteps = list(range(len(self._filenameSeries)))
                        self._currentTime = self._filenameSeries.index(filename)
                self.loadData()
                self.Modified()

    def _get_timesteps(self):
        return self._timeSteps
    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self._get_timesteps()
    def _get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self._get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

    @smproperty.stringvector(name="DataFct", information_only="1")
    def getDataFcts(self):
        return self._dataFcts
    @smproperty.stringvector(name="Datafct", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="DataFct" function="DataFct"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def setDataFcts(self, value):
        if value in self.getDataFcts():
            self._dataFct = self.getDataFcts().index(value)
            self.Modified()

    @smproperty.stringvector(name="Transform", default_values="") # , panel_visibility="never")
    @smdomain.filelist()
    @smhint.filechooser( extensions=["py"], file_description="Python script" )
    def SetTransform(self, transformPath):
        if transformPath is None:
            return
        try:
            mod = sys.modules.get(transformPath)
            if mod is None:
                mod = importlib.import_module(transformPath)
            else:
                mod = self.reload(mod)
        except ImportError:
            try:
                spec = spec_from_loader("transform", SourceFileLoader("transform", transformPath))
                mod = module_from_spec(spec)
                spec.loader.exec_module(mod)
            except FileNotFoundError:
                print("Failed to import script",transformPath)
                return
        if not hasattr(mod,"register"):
            print("Script",transformPath,"does not have a 'register' attribute - import cancelled")
            return
        transformFcts = [m.__name__ for m in mod.register]
        transformFcts[:0] = ["None"]
        transformFct  = transformFcts[0]
        self._transform = mod
        self._transformFcts = transformFcts
        self._transformFct  = transformFct
        self.Modified()

    @smproperty.stringvector(name="TransformFct", information_only="1")
    def getTransformFcts(self):
        return self._transformFcts
    @smproperty.stringvector(name="Transfct", number_of_elements="1")
    @smdomain.xml(\
        """<StringListDomain name="list">
                <RequiredProperties>
                    <Property name="TransformFct" function="TransformFct"/>
                </RequiredProperties>
            </StringListDomain>
        """)
    def setTransformFcts(self, value):
        self._transformFct = value
        self.Modified()

    @smproperty.intvector(name="Level", default_values="0")
    @smdomain.intrange(min=0, max=5)
    def SetLevel(self, level):
        self._level = level
        self.Modified()

    def loadData(self):
        ext = os.path.splitext(self._filename)[1]
        if ext == ".dgf":
            print("Still need to implement dgf reading")
            print("Which grid to use with which dimensions?")
        else:
            if self._filenameSeries is not None:
                if type(self._filenameSeries) is dict:
                    idx = self._timeSteps.index(self._currentTime)
                    with open(self._filenameSeries[str(idx)]["dumpFileName"],"rb") as f:
                        df = self.load(f)
                else:
                    with open(self._filenameSeries[self._currentTime],"rb") as f:
                        df = self.load(f)
            else:
                with open(self._filename,"rb") as f:
                    df = self.load(f)
            self._df = [d for d in df if hasattr(d,"gridView")]
            if len(self._df) > 0:
                self._gridView = self._df[0].gridView
            else:
                self._gridView = df[0]
            # make some checks:
            assert hasattr(self._gridView,"dimension"), "file read contains no valid grid view"
            assert all( [hasattr(d,"pointData") for d in self._df] ), "found a non valid grid function (no 'pointData' attribute"
            assert all( [self._gridView == d.gridView for d in self._df] ), "all grid function must be over the same gridView"
            self._dataFcts = [df.name for df in self._df]

    def RequestInformation(self, request, inInfo, outInfo):
        executive = self.GetExecutive()
        outInfo = outInfo.GetInformationObject(0)
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())

        timesteps = self._get_timesteps()
        if timesteps is not None:
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])
        return 1

    def RequestData(self, request, inInfo, outInfo):
        cTime = self._get_update_time(outInfo.GetInformationObject(0))
        if self._currentTime != cTime:
            self._currentTime = cTime
            self.loadData()
        # data
        if ( (self._transform is not None)
             and (not self._transformFct in ["","None"])
             and (not self._transformFct is None) ):
            assert self._dataFct >= 0
            gfs = getattr(self._transform,self._transformFct)\
                         (self._gridView, self._df[self._dataFct], self._df)
        else:
            gfs = self._df

        points, cells = self._gridView.tessellate(self._level)
        output = dsa.WrapDataObject(vtkUnstructuredGrid.GetData(outInfo))

        # points need to be 3d:
        if self._gridView.dimWorld == 2:
            vtk_type = vtk.VTK_TRIANGLE
            points = np.hstack([points, np.zeros((len(points), 1))])
        elif self._gridView.dimWorld == 3:
            if self._gridView.dimGrid == 2:
                vtk_type = vtk.VTK_TRIANGLE
            else:
                vtk_type = vtk.VTK_TETRA
        output.SetPoints(points)

        cell_types = np.array([], dtype=np.ubyte)
        cell_offsets = np.array([], dtype=int)
        cell_conn = np.array([], dtype=int)
        ncells, npoints = cells.shape
        cell_types = np.hstack(
                       [cell_types, np.full(ncells, vtk_type, dtype=np.ubyte)]
                     )
        offsets = len(cell_conn) + (1 + npoints) * np.arange(ncells, dtype=int)
        cell_offsets = np.hstack([cell_offsets, offsets])
        conn = np.hstack(
                   [npoints * np.ones((ncells, 1), dtype=int), cells]
               ).flatten()
        cell_conn = np.hstack([cell_conn, conn])
        output.SetCells(cell_types, cell_offsets, cell_conn)  # cell connectivities

        for df in gfs:
            array = df.pointData(self._level)
            output.PointData.append(array, df.name)  # point data
            # output.CellData.append(array, name)  # cell data
            # output.FieldData.append(array, name)  # field data

        return 1
