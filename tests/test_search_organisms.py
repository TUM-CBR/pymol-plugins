from os import path
import pytest
from tempfile import TemporaryDirectory
from tests.common import *

from cbr.sequences.visual.SearchOrganismsWidget import SearchOrganismsWidget
from cbr.extra.CbrExtraInteractiveHandler import CbrExtraInteractiveManager, run_interactive
from cbr.support.testing import PymolTestSupportObject

@pytest.fixture
def search_interactive_manager(tmp_folder: str):

    db = path.join(tmp_folder, 'test_search_organisms.sqlite')
    with run_interactive(["sequences", "interactive", f"--db-file={db}"]) as manager:
        yield manager

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
    "\t".join(["dummy", "accession", "dummy", "taxid", "name"]),
    "\t".join(["", "", "", "666",""]),
    "\t".join(["", "", "", "", "Euwallacea"])
])

class TestSearchOrganisms:
    def test_search_organisms(
        self,
        pymol_fixture: PymolTestSupportObject,
        search_organisms_widget: SearchOrganismsWidget
    ):

        pymol_fixture.run_in_ui(
            lambda: search_organisms_widget.set_search_text(dummy_organism_search_text)
        )

        from time import sleep
        sleep(5)

