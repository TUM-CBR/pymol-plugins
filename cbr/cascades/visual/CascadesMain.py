from typing import Iterable, NamedTuple, Optional, Set
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os import path
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt, pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget
import re
from typing import Any, List, Tuple

from ...core.Context import Context
from ...core.Qt.QtWidgets import with_error_handler
from ..data import CascadeStepArgs,CreateCascadeDatabaseArgs 

from .CascadesViewer import CascadesViewer
from .Ui_CascadesMain import Ui_CascadesMain

class StepEntry(NamedTuple):
    name : str
    sequences : Set[int]

email_pattern = re.compile(r"^[a-zA-Z0-9_.+-]+@[a-zA-Z0-9-]+\.[a-zA-Z0-9-.]+$")

INVALID_EMAIL_MSG = """You must provide a valid email to use this feature. This email
is provided to NCBI on every request. As this program queries the NCBI BLAST service
quite aggressively, NCBI will notify you by email if you are using the resources
too excesively.
"""

MISSING_CASCADE_ERROR = """You must define a cascade first to use this feature. Select
the FASTA file you wish to use and then select the sequences that go in every step.
"""

class CreateCascadeModel(QAbstractTableModel):

    STEP_NAME_ROW = 0
    META_ROWS = [STEP_NAME_ROW]

    SEQ_NAME_COLUMN = 0
    META_COLS = [SEQ_NAME_COLUMN]

    def __init__(
        self,
        steps : int,
        sequences : List[SeqRecord]
    ):
        super().__init__()
        self.__sequences = sequences
        self.__steps = [StepEntry(name = f"Step {i}", sequences=set()) for i in range(1, steps + 1)]

    @property
    def sequences(self):
        return self.__sequences

    def steps(self) -> List[CascadeStepArgs]:
        return [
            CascadeStepArgs(
                step_id = i,
                step_name = step.name,
                sequences = [self.__sequences[seq_ix].id for seq_ix in step.sequences]
            )
            for i,step in enumerate(self.__steps)
        ]

    def rowCount(self, parent = None) -> int:
        return len(self.META_ROWS) + len(self.__sequences)

    def columnCount(self, parent = None) -> int:
        return len(self.__steps) + len(self.META_COLS)

    def __is_step_cell(self, index: QModelIndex) -> Optional[int]:

        if index.row() == self.STEP_NAME_ROW and index.column() not in self.META_COLS:
            return index.column() - len(self.META_COLS)

    def __is_seq_name_cell(self, index: QModelIndex) -> Optional[int]:

        if index.column() == self.SEQ_NAME_COLUMN and index.row() not in self.META_ROWS:
            return index.row() - len(self.META_ROWS)

    def __is_step_checked_cell(self, index: QModelIndex) -> Optional[Tuple[int, int]]:

        if index.row() not in self.META_ROWS and index.column() not in self.META_COLS:
            return (index.column() - len(self.META_COLS), index.row() - len(self.META_ROWS))

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)
        if self.__is_step_cell(index) is not None:
            return flags | Qt.ItemIsEditable

        if self.__is_step_checked_cell(index) is not None:
            return flags | Qt.ItemIsUserCheckable

        return flags

    def set_step_count(self, steps : int):

        self.__steps = [
            self.__steps[i] if i < len(self.__steps) else StepEntry(name = f"Step {i + 1}", sequences=set())
            for i in range(0, steps)
        ]

        self.modelReset.emit()

    def __display_role_data(self, index: QModelIndex) -> Optional[str]:

        step_index = self.__is_step_cell(index)
        if step_index is not None:
            return self.__steps[step_index].name

        name_index = self.__is_seq_name_cell(index)
        if name_index is not None:
            return self.__sequences[name_index].id

    def __check_state_role(self, index: QModelIndex) -> Optional[Qt.CheckState]:

        check_state = self.__is_step_checked_cell(index)
        if check_state is None:
            return
        
        step_ix, seq_ix = check_state
        
        if seq_ix in self.__steps[step_ix].sequences:
            return Qt.Checked
        else:
            return Qt.Unchecked

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.EditRole) -> bool:

        step_ix = self.__is_step_cell(index)
        if step_ix is not None:
            self.__steps[step_ix] = self.__steps[step_ix]._replace(name = value)
            return True

        step_seq_ix = self.__is_step_checked_cell(index)
        if step_seq_ix is not None:
            step_ix, seq_ix = step_seq_ix

            if value == Qt.Checked:
                self.__steps[step_ix].sequences.add(seq_ix)
                return True
            elif value == Qt.Unchecked:
                self.__steps[step_ix].sequences.remove(seq_ix)
                return True

        return False

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.__display_role_data(index)
        if role == Qt.CheckStateRole:
            return self.__check_state_role(index)
    

class CascadesMain(QWidget):
    """This is the main widget for the Cascade BLAST application. It allows creating
    a new result set or loading an existing result set."""

    def __init__(
        self,
        context : Context
    ):
        super().__init__()
        self.__ui = Ui_CascadesMain()
        self.__context = context
        self.__ui.setupUi(self)
        self.__model = None

        self.__ui.loadExistingButton.clicked.connect(self.__on_select_db_file)
        self.__ui.selectFastaButton.clicked.connect(self.__select_fasta_file)
        self.__ui.stepsSpinBox.valueChanged.connect(self.__on_steps_changed)
        self.__ui.createButton.clicked.connect(self.__on_create_button_clicked)

    @pyqtSlot()
    def __on_steps_changed(self):

        if self.__model is None:
            return

        self.__model.set_step_count(self.__ui.stepsSpinBox.value())

    @pyqtSlot(name="__select_fasta_file")
    @with_error_handler()
    def __select_fasta_file(self):

        fasta_file,_ = QFileDialog.getOpenFileName(
            self,
            "Select a FASTA file",
            "",
            "FASTA files (*.fasta)"
        )
        
        # Open file dialog cancelled
        if fasta_file is None:
            return

        self.__set_fasta_file(SeqIO.parse(fasta_file, "fasta"))
        self.__ui.selectedFileLabel.setText(path.basename(fasta_file))

    def __set_fasta_file(self, fasta : Iterable[SeqRecord]):
        model = self.__model = CreateCascadeModel(
            self.__ui.stepsSpinBox.value(),
            list(fasta)
        )
        self.__ui.createCascadeTable.setModel(model)

    @pyqtSlot()
    def __on_select_db_file(self):
        result_file,_ = QFileDialog.getOpenFileName(
            self,
            "Select a cascade BLAST",
            "",
            "Cascade BLAST databases (*.sqlite)"
        )

        if result_file is not None:
            self.__open_db(result_file)

    def __open_db(self, db_file: str):
        self.__context.run_widget(
            lambda _: CascadesViewer(db_file)
        ).show()

    @pyqtSlot(name="__on_create_button_clicked")
    @with_error_handler()
    def __on_create_button_clicked(self):

        if self.__model is None:
            raise ValueError(
                "You must construct a cascade first!"
            )

        email = self.__ui.emailLineEdit.text()

        if email_pattern.match(email) is None:
            raise ValueError(INVALID_EMAIL_MSG)

        identity = self.__ui.identityTargetSpinBox.value() / 100

        assert 0 <= identity <= 1

        steps = self.__model.steps()

        for step in steps:
            if len(step.sequences) < 1:
                raise ValueError(
                    f"The step '{step.step_name}' has no sequences selected."
                )

        args = CreateCascadeDatabaseArgs(
            sequences = self.__model.sequences,
            target_identity = identity,
            email = email,
            steps = steps
        )

        self.__context.run_widget(
            lambda _: CascadesViewer(
                create_cascade_args=args
            )
        ).show()
