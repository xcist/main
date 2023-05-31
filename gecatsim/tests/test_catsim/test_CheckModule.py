#import pytest
import importlib.util
from gecatsim.pyfiles.CheckModules import check_module, import_module_from_spec


def test_module_not_found():
    module_spec = check_module('fake_module')
    assert module_spec is None


def test_collections_module():
    module_spec = check_module('collections')
    assert module_spec is not None
    if module_spec:
        module = import_module_from_spec(module_spec)
        print(dir(module))
        assert module is not None
        assert 'OrderedDict' in dir(module)


def test_catsim_module():
    module_spec = check_module('catsim')
    assert module_spec is not None
    if module_spec:
        module = import_module_from_spec(module_spec)
        assert module is not None
        dir_catsim = dir(module)
        print(dir_catsim)
        assert 'CatSim' in dir_catsim
        assert 'GetMu' in dir_catsim
        assert 'CFG' in dir_catsim
        assert 'rawread' in dir_catsim
        assert 'rawwrite' in dir_catsim


def test_reconstruction_module():
    module_spec = check_module('reconstruction')
    assert module_spec is not None
    if module_spec:
        module = import_module_from_spec(module_spec)
        assert module is not None
        print(dir(module))

    module_spec = check_module('fdk_equiAngle')
    assert module_spec is None
    if module_spec:
        module = import_module_from_spec(module_spec)
        print(dir(module))

    module_spec = check_module('reconstruction.pyfiles.fdk_equiAngle')
    assert module_spec is not None
    if module_spec:
        module = import_module_from_spec(module_spec)
        dir_recon_fdk_equi = dir(module)
        print(dir_recon_fdk_equi)
        assert 'fdk_equiAngle' in dir_recon_fdk_equi
        assert 'load_C_lib' in dir_recon_fdk_equi
        assert 'PathHelper' in dir_recon_fdk_equi
        assert 'mapConfigVariablesToFDK' in dir_recon_fdk_equi

    module_spec = check_module('reconstruction.pyfiles')
    assert module_spec is not None
    if module_spec:
        module = import_module_from_spec(module_spec)
        assert module is not None
        dir_recon = dir(module)
        print(dir_recon)
        assert '__package__' in dir_recon
