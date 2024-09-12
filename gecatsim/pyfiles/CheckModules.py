# from https: // www.blog.pythonlibrary.org / 2016 / 05 / 27 / python - 201 - an - intro - to - importlib /
import importlib.util


def check_module(module_name):
    """
    Checks if module can be imported without actually
    importing it
    """
    module_spec = importlib.util.find_spec(module_name)
    if module_spec is None:
        print('Module: {} not found'.format(module_name))
        return None
    else:
        print('Module: {} can be imported!'.format(module_name))
        return module_spec


def import_module_from_spec(module_spec):
    """
    Import the module via the passed in module specification
    Returns the newly imported module
    """
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module


# if __name__ == '__main__':
#     module_spec = check_module('fake_module')
#     module_spec = check_module('collections')
#     if module_spec:
#         module = import_module_from_spec(module_spec)
#         print(dir(module))
#
#     module_spec = check_module('CatSim')
#     if module_spec:
#         module = import_module_from_spec(module_spec)
#         print(dir(module))
#
#     module_spec = check_module('gecatsim.reconstruction')
#     if module_spec:
#         module = import_module_from_spec(module_spec)
#         print(dir(module))
#
#     module_spec = check_module('fdk_equiAngle')
#     if module_spec:
#         module = import_module_from_spec(module_spec)
#         print(dir(module))
#
#     module_spec = check_module('gecatsim.reconstruction.pyfiles.fdk_equiAngle')
#     if module_spec:
#         module = import_module_from_spec(module_spec)
#         print(dir(module))
#
#     module_spec = check_module('gecatsim.reconstruction.pyfiles')
#     if module_spec:
#         module = import_module_from_spec(module_spec)
#         print(dir(module))