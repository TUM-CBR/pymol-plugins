from os import path
from PyQt5.QtCore import Qt
from PyQt5.QtTest import QTest
import pytest
from tests.common import *

from cbr.sequences.visual.SearchOrganismsWidget import SearchOrganismsWidget
from cbr.extra.CbrExtraInteractiveHandler import CbrExtraInteractiveManager, run_interactive
from cbr.support.testing import PymolTestSupportObject

@pytest.fixture
def search_interactive_manager(tmp_folder: str, pymol_fixture: PymolTestSupportObject):

    db = path.join(tmp_folder, 'test_search_organisms.sqlite')

    def create():
        return run_interactive(["sequences", "interactive", f"--db-file={db}"])

    manager = pymol_fixture.run_in_ui(create)
    yield manager
    pymol_fixture.run_in_ui(lambda: manager.stop())

@pytest.fixture
def search_organisms_widget(
    pymol_fixture: PymolTestSupportObject,
    search_interactive_manager: CbrExtraInteractiveManager
):
    def create_widget():
        widget = SearchOrganismsWidget(search_interactive_manager)
        widget.show()

        # The window must be created before this function returns
        import time
        time.sleep(1)
        return widget

    yield pymol_fixture.run_in_ui(create_widget)

dummy_organism_search_text = "\n".join([
    "\t".join(["dummy1", "accession", "dummy2", "taxid", "name"]),
    "\t".join(["attribute1", "", "attributex", "666",""]),
    "\t".join(["attribute2", "", "attributey", "", "Euwallacea"])
])

class TestSearchOrganisms:

    def test_1_empty_search_organisms(
        self,
        search_organisms_widget: SearchOrganismsWidget,
        pyqt_test_helpers: PyQtTestHelpers
    ):
        QTest.mouseClick(search_organisms_widget.search_button, Qt.LeftButton)
        with pytest.raises(Exception, match="ValueError"):
            pyqt_test_helpers.wait_in_ui([search_organisms_widget], lambda: False)

    def test_2_search_organisms(
        self,
        search_organisms_widget: SearchOrganismsWidget,
        pyqt_test_helpers: PyQtTestHelpers
    ):

        pyqt_test_helpers.run_in_ui(
            lambda: search_organisms_widget.set_search_text(dummy_organism_search_text)
        )

        def select_columns():
            search_organisms_widget.name_combo.setCurrentText("name")
            assert search_organisms_widget.name_combo.currentText() == "name"

            search_organisms_widget.taxid_combo.setCurrentText("taxid")
            assert search_organisms_widget.taxid_combo.currentText() == "taxid"

        pyqt_test_helpers.run_in_ui(select_columns)

        import time
        time.sleep(1)
        QTest.mouseClick(search_organisms_widget.search_button, Qt.LeftButton)

        def has_results():
            model = search_organisms_widget.results_table.model()
            for row in range(model.rowCount()):
                for col in range(model.columnCount()):
                    item = model.data(model.index(row, col))
                    if item != "":
                        return True

            return False

        pyqt_test_helpers.wait_in_ui([search_organisms_widget], has_results, timeout=5)

        results_model = search_organisms_widget.results_table.model()

        for row in range(results_model.rowCount()):
            for col in range(results_model.columnCount()):
                item = results_model.data(results_model.index(row, col))
                print(item)
