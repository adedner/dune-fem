from argparse import ArgumentParser

from dune.source.cplusplus import NameSpace, SourceWriter

from .formfiles import compileUFLFile

def compile(inputFileName, outputFileName=None, namespace=None, tempVars=True):
    code = compileUFLFile(inputFileName, tempVars=tempVars)
    if namespace is not None:
        for ns in reversed(namespace.split('::')):
            code = NameSpace(ns, code=code)
    if outputFileName is None:
        if inputFileName[-4:] != ".ufl":
            outputFileName = inputFileName + ".hh"
        else:
            outputFileName = inputFileName[:-4] + ".hh"
    SourceWriter(outputFileName).emit(code)


def version(package_name):
    try:
        import pip
        from pip._internal.utils.misc import get_installed_distributions
        packages = get_installed_distributions()
    except ImportError:
        'pip not found, skipping version check'
        return "?.?"
    except:
        'using pip version < 10'
        packages = pip.get_installed_distributions()
    for package in packages:
        if package.project_name == package_name and package.has_version():
            return str(package.version)
    return "?.?"


def main():
    description = 'dune-fempy conservationlaw model compiler'

    parser = ArgumentParser(description=description)
    parser.add_argument('input', help='name of input .ufl file')
    parser.add_argument('-o', '--output', help='name of output .hh file')
    parser.add_argument('-n', '--namespace', help='C++ namespace for models in')
    parser.add_argument('--no-tempvars', dest='tempVars', action='store_false', help='do not generate temporary variables')
    parser.add_argument('-v', '--version', action='version', version=description + " " + version("dune.fem"))
    parser.add_argument('--debug', action='store_true', help='do not catch exceptions')
    try:
        args = parser.parse_args()
    except Exception as e:
        print(e)
        return 1

    try:
        compile(args.input, outputFileName=args.output, namespace=args.namespace, tempVars=args.tempVars)
    except Exception as e:
        if args.debug:
            raise(e)
        else:
            print(e)
            return 1

    return 0
