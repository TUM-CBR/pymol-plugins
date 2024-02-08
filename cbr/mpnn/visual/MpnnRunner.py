from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QShowEvent
from PyQt5.QtWidgets import QWidget
from pymol import cmd
from typing import Dict, Iterable, List, NamedTuple, Optional, Set, Tuple

from ...core.Context import Context, KnownExecutables
from ...core.Qt.QtWidgets import show_exception, with_error_handler
from ...core.Qt.visual.NamedTupleEditor import MetaFieldOverrides, namedtuple_eidtor
from ..data import MpnnEditSpace, MpnnSpec
from .MpnnViewer import MpnnViewer
from .Ui_MpnnRunner import Ui_MpnnRunner

def with_selection(
    current: MpnnSpec,
    selections : Iterable[str]
) -> MpnnSpec:

    results : Dict[Tuple[str, str], MpnnEditSpace] = {}

    def add_item(model: str, chain: str, resv: int) -> None:
        key = (model, chain)
        current = results.get(key)

        if current is None:
            results[key] = current = MpnnEditSpace(
                model=model,
                chain=chain,
                residues=set()
            )

        current.residues.add(resv)

    for selection in selections:
        # We only want the residue position of the selection
        # for that reason we use polymer.protein
        cmd.iterate(
            f"{selection} & polymer.protein",
            'add_item(model, chain, resv)',
            space={'add_item': add_item}
        )

    return current._replace(
        edit_spaces=list(results.values())
    )

class MpnnSelection(NamedTuple):
    name : str
    include: bool

class MpnnRunner(QWidget):

    def __init__(
        self,
        context: Context
    ):
        super().__init__()
        self.__context = context
        self.__ui = Ui_MpnnRunner()
        self.__ui.setupUi(self)
        self.__ui.runButton.clicked.connect(self.__on_run_clicked)
        self.__ui.refreshButton.clicked.connect(self.__on_refresh_button_clicked)
        self.__selections_model = namedtuple_eidtor(
            self.__ui.selectionsTable,
            tuple_type=MpnnSelection,
            tuple_field_overrides = {
                'name': MetaFieldOverrides(readonly=True)
            }
        )

    def __is_valid(self, spec: MpnnSpec) -> Optional[str]:
        canary = False

        for space in spec.edit_spaces:
            canary = True
            if len(space.residues) == 0:
                return "Your selection does not contain any residues"
            
        if not canary:
            return "You must select residues to use this feature"
        
    @pyqtSlot()
    def __on_refresh_button_clicked(self):
        self.__refresh_selections()
        
    def __refresh_selections(self):
        selections : Set[str] = set(cmd.get_names('selections'))
        removed : List[MpnnSelection] = []
        model = self.__selections_model

        current_items = list(model.iterate_values())

        for item in current_items:
            
            assert item is not None, "All items must have a value"

            # Item is not in selections, meaning that the selection
            # no longer exists.
            if item.name not in selections:
                removed.append(item)

            # REmove item from selections as it is already
            # being displayed by the table
            else:
                selections.remove(item.name)

        model.remove(*removed)
        model.append(*
            (
                MpnnSelection(name=selection, include=False)
                for selection in selections
            )
        )

    def __ensure_components(self):
        
        try:
            # Enusre the ProteinMPNN executable is available
            self.__context.get_executable(self, KnownExecutables.ProteinMPNN)
        except Exception as e:
            show_exception(self, e)
            self.close()
        
    def showEvent(self, a0: Optional[QShowEvent]):

        if not self.isVisible():
            return

        self.__ensure_components()
        self.__refresh_selections()

    @pyqtSlot(name="__on_run_clicked")
    @with_error_handler()
    def __on_run_clicked(self):

        spec = MpnnSpec(
            num_seqs=self.__ui.numSequencesSpinBox.value(),
            edit_spaces=[]
        )

        spec = with_selection(
            spec,
            (
                item.name
                for item in self.__selections_model.iterate_values()
                    if item is not None and item.include
            )
        )
        validate = self.__is_valid(spec)

        if validate is not None:
            raise Exception(validate)

        mpnn = self.__context.get_executable(self, KnownExecutables.ProteinMPNN)
        self.__context.run_widget(
            lambda ctx: MpnnViewer(ctx, mpnn, spec)
        ).show()