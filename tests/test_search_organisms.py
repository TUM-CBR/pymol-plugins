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

        return widget

    yield pymol_fixture.run_in_ui(create_widget)

dummy_organism_search_text = "\n".join([
    "\t".join(["dummy1", "accession", "dummy2", "taxid", "name"]),
    "\t".join(["attribute1", "", "attributex", "666",""]),
    "\t".join(["attribute2", "", "attributey", "", "Oxyura"])
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
            return search_organisms_widget.busy_progress.isHidden()

        pyqt_test_helpers.wait_in_ui([search_organisms_widget], has_results, timeout=100)

        results_model = search_organisms_widget.results_table.model()

        expected_accessions = ["GCA_019458465.1", "GCF_003574155.1", "GCA_011077185.1", "GCF_011077185.1"]
        expected_attributes = ["attribute1", "attribute2", "attributex", "attributey"]

        def check_results():
            for row in range(results_model.rowCount()):
                for col in range(results_model.columnCount()):
                    text = results_model.data(results_model.index(row, col))
                    if text is None:
                        continue

                    if text in expected_accessions:
                        expected_accessions.pop(expected_accessions.index(text))
                    
                    if text in expected_attributes:
                        expected_attributes.pop(expected_attributes.index(text))

        pyqt_test_helpers.run_in_ui(check_results)
        assert len(expected_accessions) == 0, f"Missing accessions: {expected_accessions}"
        assert len(expected_attributes) == 0, f"Missing attributes: {expected_attributes}"

    @pytest.mark.skip("This test can be used to run the widget in isolation.")
    def test_3_demo(
        self,
        search_organisms_widget: SearchOrganismsWidget,
        pyqt_test_helpers: PyQtTestHelpers
    ):
        input()
