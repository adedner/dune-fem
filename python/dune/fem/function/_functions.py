from __future__ import absolute_import, division, print_function, unicode_literals

import sys, os
import logging
logger = logging.getLogger(__name__)

import dune.grid
import dune.fem.space
import dune.models.localfunction

import dune.common.checkconfiguration as checkconfiguration
from dune.common.hashit import hashIt

def registerGridFunctions(gridView):
    from dune.generator import builder
    typeName = gridView._typeName
    moduleName = "femgridfunctions_" + hashIt(typeName)

    includes = ["dune/fempy/py/grid/gridpart.hh", "dune/fempy/py/grid/function.hh"] + gridView._includes

    source = "#include <config.h>\n\n"
    source += "".join(["#include <" + i + ">\n" for i in includes])
    source += "\n"
    source += "PYBIND11_MODULE( " + moduleName + ", module )\n"
    source += "{\n"
    source += "  typedef Dune::FemPy::GridPart< " + gridView._typeName + "> GridPart;\n"
    source += '  module.def( "globalGridFunction", Dune::FemPy::defGlobalGridFunction< GridPart >( module, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ));\n'
    source += '  module.def( "localGridFunction", Dune::FemPy::defLocalGridFunction< GridPart > ( module, "LocalGridFunction",  std::make_integer_sequence< int, 11 >() ));\n'
    source += "}\n"

    return builder.load(moduleName, source, "gridfunctions")

def globalFunction(gridView, name, order, value):
    # assert False
    module = registerGridFunctions(gridView)
    return module.globalGridFunction(gridView,name,order,value).as_ufl()


def localFunction(gridView, name, order, value):
    # assert False
    module = registerGridFunctions(gridView)
    return module.localGridFunction(gridView,name,order,value).as_ufl()


def levelFunction(gridView,name="levels"):
    @dune.grid.gridFunction(gridView,name=name)
    def levelFunction(e,x):
        return [e.level]
    return levelFunction


def partitionFunction(gridView,name="rank"):
    @dune.grid.GridFunction(gridView,name=name)
    class Partition(object):
        def __init__(self,rank):
            self.rank = rank
        def __call__(self,en,x):
            return [self.rank]
    return Partition(gridView.comm.rank) # localFunction(gridView, name, 0, Partition(gridView.comm.rank))


def cppFunction(gridView, name, order, code, *args, **kwargs):
    raise NotImplementedError("cpp function is not working at the moment")
    return dune.models.localfunction.generatedFunction(gridView, name, order, code, *args, **kwargs)


def uflFunction(gridView, name, order, ufl, virtualize=True, scalar=False, *args, **kwargs):
    Func = dune.models.localfunction.UFLFunction(gridView, name, order,
            ufl, renumbering=None,
            virtualize=virtualize, *args, **kwargs)
    if Func is None:
        raise AttributeError("could not generate ufl grid function from expression "+str(ufl))
    func = Func(gridView,name,order,*args,**kwargs)
    if not hasattr(func,"scalar"):
        func.scalar = scalar
    return func.as_ufl() if func is not None else None

def discreteFunction(space, name, expr=None, dofVector=None):
    """create a discrete function

    Args:
        space: discrete function space
        name:  name of the discrete function
        expr:  analytical expression to interpolate

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    if dofVector is None:
        df = space.DiscreteFunction(space,name)
    else:
        df = space.DiscreteFunction(name,space,dofVector)

    if expr is None and dofVector is None:
        df.clear()
    elif expr is not None:
        df.interpolate(expr)
    df.scalar = space.scalar
    return df.as_ufl()

def tupleDiscreteFunction(*spaces, **kwargs):
    # from dune.fem.discretefunction import module, addAttr
    name = kwargs.get("name", "")
    try:
        tupleSpace = spaces[0]
        spaces = spaces[0].components
    except AttributeError:
        tupleSpace = dune.fem.space.tuple(*spaces)
    DiscreteFunction = tupleSpace.DiscreteFunction
    df = DiscreteFunction(tupleSpace,name)
    compNames = kwargs.get("components", None)
    if not compNames is None:
        components = df.components
        assert len(compNames) == len(components)
        for c, n in zip(components, compNames):
            setattr(df,n,c)
    return df.as_ufl()
