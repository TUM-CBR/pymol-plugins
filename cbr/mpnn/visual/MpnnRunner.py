from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QShowEvent
from PyQt5.QtWidgets import QWidget
from pymol import cmd
from typing import Dict, Iterable, List, NamedTuple, Optional, Set, Tuple


from ...core.Context import Context, KnownExecutables
from ...core.pymol.visual.PymolChainSetEditor import PymolChainSetEditorModel
from ...core.Qt.QtWidgets import show_exception, with_error_handler
from ...core.Qt.visual.NamedTupleEditor import FieldOrientation, MetaFieldOverrides, namedtuple_eidtor
from ...core.sequence import assert_residue
from ..data import MpnnEditSpace, MpnnSpec, models, TiedPositionsSpec
from .MpnnAdvancedOptionsDialog import MpnnAdvancedOptionsDialog
from .MpnnViewer import MpnnViewer
from .Ui_MpnnRunner import Ui_MpnnRunner

class MpnnSelection(NamedTuple):
    """
    This class specifies a PyMol selection which represents a section of
    the protien that will be designed.

    Attribtues:
        name                The name of the pymol selection containing the residues to be designed.
        include             Whether that section will be designed.
        exculded_residues   What residues to exclude when designing that section. Sequence of single letters.
    """

    name : str
    include: bool
    excluded_residues: str = ""

    def get_exclude_list(self) -> Optional[Set[str]]:
        excluded = self.excluded_residues.strip()

        result = set(
            assert_residue(letter)
            for letter in excluded
                if letter.isalpha()
        )

        if len(result) == 0:
            return None
        
        return result

def with_selection(
    current: MpnnSpec,
    selections : Iterable[MpnnSelection]
) -> MpnnSpec:

    results : Dict[Tuple[str, str], MpnnEditSpace] = {}

    for selection in selections:

        excluded = selection.get_exclude_list()
        def add_item(model: str, chain: str, resv: int) -> None:
            key = (model, chain)
            current = results.get(key)

            if current is None:
                results[key] = current = MpnnEditSpace(
                    model=model,
                    chain=chain,
                    residues=set(),
                    excluded={}
                )

            current.residues.add(resv)
            if excluded is not None:
                current.excluded[resv] = excluded

    
        # We only want the residue position of the selection
        # for that reason we use polymer.protein
        cmd.iterate(
            f"{selection.name} & polymer.protein",
            'add_item(model, chain, resv)',
            space={'add_item': add_item}
        )

    return current._replace(
        edit_spaces=list(results.values())
    )

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
            },
            fields_orientation=FieldOrientation.Horizontal
        )
        self.__advanced_settings = MpnnAdvancedOptionsDialog()
        self.__ui.advancedOptionsButton.clicked.connect(self.__on_advanced_options_clicked)
        self.__chain_set_model = PymolChainSetEditorModel(exclusive=True)
        self.__ui.tiedPositionsTable.setModel(self.__chain_set_model)
        self.__ui.modelComboBox.addItems(models.keys())

    @pyqtSlot()
    def __on_advanced_options_clicked(self):
        self.__advanced_settings.exec_()

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

        self.__chain_set_model.refresh()
        self.__ui.tiedPositionsTable.resizeColumnsToContents()

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
            edit_spaces=[],
            mpnn_args=self.__advanced_settings.value(),
            tied_positions=TiedPositionsSpec(tied_chains=self.__chain_set_model.values()),
            mpnn_model=models[self.__ui.modelComboBox.currentText()]
        )

        spec = with_selection(
            spec,
            (
                item
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