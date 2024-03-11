import json
from os import path
from pymol import cmd
from PyQt5.QtCore import QObject
from PyQt5.QtWidgets import QWidget
from typing import Any, Dict, NamedTuple, Optional

from ...core.Context import Context
from ...core.pymol.structure import StructureSelection
from ...core.Qt.QtWidgets import show_exception, with_error_handler
from ...core.uRx import dsl as rx
from ...core.uRx import qt as qtRx
from ...core.uRx import util as rxUtil
from ...core.visual import as_structure_selector
from ...extra.CbrExtraInteractiveHandler import CbrExtraInteractiveHandler, CbrExtraInteractiveManager, CbrProcessExit, MessageQueingPolicy, run_interactive
from ..data import CavitiesInteractiveInput, CavitiesInteractiveOutput, FindCavitiesArgs
from .Ui_CavityFinder import Ui_CavityFinder

K_CAVITIES = "cavities"
K_POINTS = "points"
K_RADII = "radii"
K_ID = "id"

def structure_id(structure: StructureSelection) -> str:
    return structure.show()

def create_structure_json(
    structrure: StructureSelection,
    tmp_dir: str,
    state: int = 1
):
    coords = []
    vdw_radii = []

    def apply(x: float, y: float, z: float, vdw: float):
        coords.append([x,y,z])
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
        K_ID: structure_id(structrure),
        K_POINTS: coords,
        K_RADII: vdw_radii
    }

    json_file_path = path.join(tmp_dir, f"{id(json_object)}.json")

    with open(json_file_path, 'w') as json_stream:
        json.dump(
            [json_object],
            json_stream
        )

    return json_file_path

class CavityFinderInstance(NamedTuple):
    interactive_manager: CbrExtraInteractiveManager
    structure: StructureSelection
    cavity_handler: CbrExtraInteractiveHandler[CavitiesInteractiveOutput, CavitiesInteractiveInput]

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
        cavity_handler = interactive.message_handler(
            CavitiesInteractiveOutput.parse,
            CavitiesInteractiveInput.serialize,
            MessageQueingPolicy.HANDLE_LATEST
        )

        return CavityFinderInstance(interactive, structure, cavity_handler)

class CavityFinder(QWidget):

    def __init__(self, context: Context) -> None:
        super().__init__()

        self.__ui = Ui_CavityFinder()
        self.__ui.setupUi(self)

        self.__structure_selector = as_structure_selector(
            self.__ui.structureCombo,
            self.__ui.refreshButton
        )

        self.__click_os = qtRx.observer_from_signal0(self, self.__ui.findButton.clicked, hot=False)
        self.__busy_subscription = rxUtil.UpdatableSubscription()
        self.__result_subscription = rxUtil.UpdatableSubscription()
        self.__status_subscription = rxUtil.UpdatableSubscription()
        rx.observe(self.__click_os).for_each(
            lambda _: self.__on_find_button_clicked()
        )
        self.__cavity_process: Optional[CavityFinderInstance] = None
        self.__working_dir = context.create_temporary_directory()
        self.__ui.busyProgress.hide()

    def __get_or_create_cavity_process(self):
        structure = self.__structure_selector.currentSelection
        cavity_process = self.__cavity_process

        if structure is None:
            raise ValueError("No structure has been selected!")

        if cavity_process is None or cavity_process.structure != structure:
            cavity_process.stop() if cavity_process is not None else None
            cavity_process = CavityFinderInstance.start(self, structure, self.__working_dir)

            # True means the app is busy, false means not busy
            self.__busy_subscription.update(
                rx.observe(self.__click_os).map(lambda _: True) \
                    .merge(
                        cavity_process.cavity_handler.observe_values().map(
                            lambda _: False
                        )
                    ).for_each(
                        self.__on_busy_changed,
                        on_error=self.__on_error
                    )
            )

            self.__result_subscription.update(
                cavity_process.cavity_handler.observe_values() \
                    .for_each(
                        self.__on_result,

                        # Errors are handled by the previous subscription
                        on_error=lambda _: None
                    )
            )

            self.__status_subscription.update(
                cavity_process.interactive_manager.observe_status() \
                    .for_each(
                        self.__on_quit,
                        on_error=self.__on_error
                    )
            )

        return cavity_process
    
    def __on_quit(self, _result: CbrProcessExit):
        pass
    
    def __on_result(self, result: CavitiesInteractiveOutput):
        print(result)

    def __on_error(self, exn: Exception):
        show_exception(self, exn)
        self.__on_busy_changed(False)

    def __on_busy_changed(self, is_busy: bool) -> None:
        self.__ui.findButton.setEnabled(not is_busy)
        self.__ui.busyProgress.setVisible(is_busy)

    def __get_cavity_args(self, structure_id: str) -> FindCavitiesArgs:
        min_volume = self.__ui.minSizeSelector.value()
        max_volume = self.__ui.maxSizeSelector.value()

        return FindCavitiesArgs(
            points_id=structure_id,
            min_volume=min_volume,
            max_volume=max_volume
        )

    @with_error_handler()
    def __on_find_button_clicked(self):

        processor = self.__get_or_create_cavity_process()
        args = self.__get_cavity_args(structure_id(processor.structure))
        processor.cavity_handler.send_message(
            CavitiesInteractiveInput(
                find_cavities=[args]
            )
        )