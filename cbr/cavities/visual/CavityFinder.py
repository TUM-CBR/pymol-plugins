import json
from os import path
from pymol import cmd
from PyQt5.QtCore import pyqtSlot, QObject, Qt
from PyQt5.QtWidgets import QWidget
import tempfile
from typing import Any, Dict, NamedTuple, Optional

from ...core.Context import Context
from ...core.pymol.structure import StructureSelection
from ...core.Qt.QtWidgets import with_error_handler
from ...core.visual import as_structure_selector
from ...extra.CbrExtraInteractiveHandler import CbrExtraInteractiveHandler, CbrExtraInteractiveManager, run_interactive
from .Ui_CavityFinder import Ui_CavityFinder

K_CAVITIES = "cavities"
K_POINTS = "points"
K_RADII = "radii"

def create_structure_json(
    structrure: StructureSelection,
    tmp_dir: str,
    state: int = 1
):
    coords = []
    vdw_radii = []

    def apply(x: float, y: float, z: float, vdw: float):
        coords.append([coords])
        vdw_radii.append(vdw)

    cmd.iterate_state(
        state,
        structrure.selection,
        'apply(x,y,z,vdw)',
        space={'apply': apply}
    )

    # This is the type that describes this json object
    # https://github.com/TUM-CBR/cbr-tools-extra/blob/d2b1b573e0216695ced93fd3c434966de6f030e1/cbrextra/cavities/data.py#L35
    json_object: Dict[Any, Any] = {
        K_CAVITIES: {
            structrure.show(): {
                K_POINTS: coords,
                K_RADII: vdw_radii
            }
        }
    }

    json_file_path = path.join(tmp_dir, f"{id(json_object)}.json")

    with open(json_file_path, 'w') as json_stream:
        json.dump(
            json_object,
            json_stream
        )

    return json_file_path

class CavityFinderInstance(NamedTuple):
    interactive_manager: CbrExtraInteractiveManager
    structure: StructureSelection

    def stop(self):
        pass

    @classmethod
    def start(
        cls,
        parent: QObject,
        structure: StructureSelection,
        working_folder: str
    ) -> 'CavityFinderInstance':
        
        json_file = create_structure_json(structure, working_folder)
        interactive = run_interactive(
            [
                "cavities",
                "interactive",
                "--input-points",
                json_file
            ],
            parent=parent
        )

        return CavityFinderInstance(interactive, structure)

class CavityFinder(QWidget):

    def __init__(self, context: Context) -> None:
        super().__init__()

        self.__ui = Ui_CavityFinder()
        self.__ui.setupUi(self)

        self.__structure_selector = as_structure_selector(
            self.__ui.structureCombo,
            self.__ui.refreshButton
        )

        self.__ui.findButton.clicked.connect(self.__on_find_button_clicked)
        self.__cavity_process: Optional[CavityFinderInstance] = None
        self.__working_dir = context.create_temporary_directory()

    def get_or_create_cavity_process(self):
        structure = self.__structure_selector.currentSelection
        cavity_process = self.__cavity_process

        if structure is None:
            raise ValueError("No structure has been selected!")

        if cavity_process is None or cavity_process.structure != structure:
            cavity_process.stop() if cavity_process is not None else None
            cavity_process = CavityFinderInstance.start(self, structure, self.__working_dir)

        return cavity_process

    @pyqtSlot()
    @with_error_handler()
    def __on_find_button_clicked(self):

        pass
