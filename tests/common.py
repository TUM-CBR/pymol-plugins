from PyQt5.QtCore import QObject
from PyQt5.QtWidgets import QWidget
import pytest
import time
from typing import Callable, Dict, List, Optional

@pytest.fixture(scope="session")
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

@pytest.fixture(scope="function")
def tmp_folder():
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmpdirname:
        yield tmpdirname

class UiErrorsException(Exception):
    def __init__(self, errors: List[Exception]):
        super().__init__(
            "Many errors occurred: \n" + "\n".join(
                [
                    f"{error}"
                    for error in errors
                ]
            )
        )
        self.errors = errors

class PyQtTestHelpers(QObject):

    def __init__(
        self,
        pymol_support
    ):
        super().__init__()
        self.__pymol_support = pymol_support
        self.__exns: Dict[int, List[Exception]] = {}
    
    def show_error_patched(
        self,
        parent: Optional[QWidget],
        title: str,
        description: str = '',
        window_title: str = 'Error'
    ):

        uid = parent.__hash__() if parent is not None else 0
        exns = self.__exns.get(uid)

        if exns is None:
            exns = self.__exns[uid] = []

        exns.append(Exception(title, description))

    def __collet_errors(self, parent: Optional[QWidget]) -> List[Exception]:
        uid = parent.__hash__() if parent is not None else 0
        return self.__exns.pop(uid, [])

    def __collect_all_errors(self, parents: List[QWidget]) -> List[Exception]:
        return [
            error
            for parent in parents
            for error in self.__collet_errors(parent)
        ]

    def __raise_if_errors(self, parents: List[QWidget]):
        errors = self.__collect_all_errors(parents)
        if len(errors) == 1:
            raise errors[0]
        elif len(errors) > 1:
            raise UiErrorsException(errors)

    def wait_in_ui(
        self,
        widgets: List[QWidget],
        priedicate: Callable[[], bool],
        timeout: int = 10
    ):

        start = time.time()

        while time.time() - start < timeout:
            if self.__pymol_support.run_in_ui(priedicate):
                self.__raise_if_errors(widgets)
                return

            self.__raise_if_errors(widgets)

            time.sleep(0.1)

        raise TimeoutError("The predicate did not become ture in the alloted time")

@pytest.fixture(scope="function")
def pyqt_test_helpers(monkeypatch, pymol_fixture):
    helpers = PyQtTestHelpers(pymol_fixture)
    monkeypatch.setattr("pmg_tk.startup.cbr.core.Qt.QtWidgets.show_error", helpers.show_error_patched)
    monkeypatch.setattr("cbr.core.Qt.QtWidgets.show_error", helpers.show_error_patched)
    return helpers

