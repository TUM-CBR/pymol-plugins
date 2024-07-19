import pytest

@pytest.fixture
def pymol_fixture():
    import pymol
    import time
    pymol.finish_launching()

    pymol_support = None

    while pymol_support is None:
        try:
            from pmg_tk.startup.cbr.support.testing import PymolTestSupportObject
            pymol_support = PymolTestSupportObject()
        except ImportError:
            time.sleep(0.1)

    with pymol_support:
        yield pymol_support
