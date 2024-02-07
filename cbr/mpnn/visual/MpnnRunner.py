from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QWidget
from pymol import cmd
from typing import Dict, Optional, Tuple

from ...core.Context import Context, KnownExecutables
from ...core.Qt.QtWidgets import with_error_handler
from ..data import MpnnEditSpace, MpnnSpec
from .MpnnViewer import MpnnViewer
from .Ui_MpnnRunner import Ui_MpnnRunner

def with_selection(
    current: MpnnSpec,
    selection : str = 'sele'
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

    cmd.iterate(
        selection,
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
        self.__init_widget()

    def __is_valid(self, spec: MpnnSpec) -> Optional[str]:
        canary = False

        for space in spec.edit_spaces:
            canary = True
            if len(space.residues) == 0:
                return "Your selection does not contain any residues"
            
        if not canary:
            return "You must select residues to use this feature"
        
    def __init_widget(self):

        # Enusre the ProteinMPNN executable is available
        self.__context.get_executable(self, KnownExecutables.ProteinMPNN)

    @pyqtSlot(name="__on_run_clicked")
    @with_error_handler()
    def __on_run_clicked(self):

        spec = MpnnSpec(
            num_seqs=self.__ui.numSequencesSpinBox.value(),
            edit_spaces=[]
        )

        spec = with_selection(spec)
        validate = self.__is_valid(spec)

        if validate is not None:
            raise Exception(validate)

        mpnn = self.__context.get_executable(self, KnownExecutables.ProteinMPNN)
        self.__context.run_widget(
            lambda ctx: MpnnViewer(ctx, mpnn, spec)
        )