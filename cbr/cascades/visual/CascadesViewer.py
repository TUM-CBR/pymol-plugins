from Bio import SeqIO
from os import path
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, QProcess, Qt, pyqtSlot
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QFileDialog, QWidget
import shutil
import tempfile
from typing import Any, Callable, Dict, Iterable, Optional, Tuple, TypeVar

from ...core import color
from ...core.Qt.QtWidgets import show_error, show_info, with_error_handler
from ..data import *
from ..operations import query_organisms, create_cascade

from .CascadeFilterWidget import CascadeFilterWidget
from .CascadeFilterWidget import CascadeFilterWidget, CascadeFilterOperator
from .ProgressWidget import ProgressWidget
from .Ui_CascadesViewer import Ui_CascadesViewer

TOut = TypeVar('TOut')

identity_color_spread = color.color_spread(
    (100,0,0),
    (0,100,0)
)

white = QColor(255,255,255)

class OrganismsTableModel(QAbstractTableModel):

    POLICY_ROW_IX = 0
    META_ROWS = [POLICY_ROW_IX]

    ORGANISM_NAME_COL_IX = 0
    META_COLS = [ORGANISM_NAME_COL_IX]
    META_COLS_HEADERS = ["Organism"]

    def __init__(
        self,
        results : QueryCascadeResult
    ) -> None:
        super().__init__()
        self.__organism_results = results
        self.__step_policies = {}

    def __get_policy(self, step_id : int):
        step = self.steps[step_id]
        return self.__step_policies.get(step.step_id) or QueryStepPolicy.any

    @property
    def __organisms(self):
        return self.__organism_results.organisms

    def enumerate_organisms_rows(self) -> Iterable[Tuple[int, QueryCascadeResultOrganismEntry]]:
        organisms = self.__organisms
        offset = len(self.META_ROWS)
        return zip(
            range(offset, offset + len(organisms)),
            organisms
        )

    @property
    def steps(self) -> List[QueryCascadeResultStepEntry]:
        if len(self.__organisms) > 0:
            return self.__organisms[0].steps
        else:
            return []

    def __get_headers(self):
        return self.META_COLS_HEADERS + [step.step_name for step in self.steps]

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            headers = self.__get_headers()
            if 0 <= section < len(headers):
                return headers[section]
        return super().headerData(section, orientation, role)

    def columnCount(self, parent = None) -> int:
        return len(self.META_COLS) + len(self.steps)

    def rowCount(self, parent = None) -> int:
        return len(self.__organisms) + len(self.META_ROWS)

    def __get_organism(self, index: QModelIndex) -> QueryCascadeResultOrganismEntry:
        ix = index.row() - len(self.META_ROWS)
        return self.__organisms[ix]

    def __get_organism_record(
        self,
        col_ix : int,
        entry : QueryCascadeResultOrganismEntry,
        with_organism : Callable[[OrganismResultEntry], TOut],
        with_step : Callable[[QueryCascadeResultStepEntry], TOut]
    ) -> TOut:
        if col_ix == 0:
            return with_organism(entry.organism)
        else:
            step = entry.steps[col_ix - 1]
            return with_step(step)

    def __render_organism(
        self,
        col_ix : int,
        entry : QueryCascadeResultOrganismEntry
    ) -> str:
        return self.__get_organism_record(
            col_ix,
            entry,
            with_organism=lambda o: f"{o.name} (taxid:{o.tax_id})",
            with_step=lambda s: f"{round(s.identity*100,2)}%"
        )

    def __display_role_data(self, index : QModelIndex) -> Any:
        row_ix = index.row()
        col_ix = index.column()

        if row_ix == self.POLICY_ROW_IX:

            meta_cols_count = len(self.META_COLS)
            if col_ix < meta_cols_count:
                return ""
            
            return self.__get_policy(col_ix - meta_cols_count)

        return self.__render_organism(col_ix, self.__get_organism(index))

    def __color_role_data(self, index : QModelIndex) -> Any:

        if index.row() < len(self.META_ROWS) or index.column() < len(self.META_COLS):
            return None

        organism = self.__get_organism(index)
        return self.__get_organism_record(
            index.column(),
            organism,
            with_organism=lambda _: None,
            with_step=lambda step: QColor(*identity_color_spread.get_color(step.identity))
        )

    def flags(self, index : QModelIndex) -> Qt.ItemFlags:
        flags = super().flags(index)

        if index.row() == self.POLICY_ROW_IX and index.column() > 1:
            return flags | Qt.ItemIsEditable

        return flags

    def data(self, index: QModelIndex, role: int = Qt.DisplayRole) -> Any:

        if not index.isValid():
            return None

        if role == Qt.DisplayRole:
            return self.__display_role_data(index)
        elif role == Qt.BackgroundColorRole:
            return self.__color_role_data(index)
        elif role == Qt.TextColorRole:
            return self.__color_role_data(index) and white

        return None

class CascadesViewer(QWidget):

    def __init__(
        self,
        db_path : Optional[str] = None,
        create_cascade_args : Optional[CreateCascadeDatabaseArgs] = None
    ):
        super().__init__()

        self.__ui = Ui_CascadesViewer()
        self.__ui.setupUi(self)
        self.__model = OrganismsTableModel(QueryCascadeResult(organisms=[]))
        self.__policy_combos : Dict[int, CascadeFilterWidget] = {}
        self.__cascade_process = None
        self.__progress_widget = ProgressWidget()

        self.layout().replaceWidget(
            self.__ui.progressWidget,
            self.__progress_widget
        )

        self.__action_widgets : List[QWidget] = [
            self.__ui.queryButton,
            self.__ui.saveResultsButton
        ]

        self.__toggle_working_widgets(False)
        self.__ui.queryButton.clicked.connect(self.__filter_results)
        self.__ui.saveResultsButton.clicked.connect(self.__save_database)

        if db_path is None:
            self.__tmp_dir = tempfile.TemporaryDirectory()
            self.__db_path = path.join(self.__tmp_dir.name, "cascades.sqlite")
        else:
            self.__tmp_dir = None
            self.__db_path = db_path

        if create_cascade_args is not None:
            self.__create_cascade(create_cascade_args)
        else:
            self.__create_initial_table()

    @pyqtSlot()
    def __save_database(self):
        
        result, _ = QFileDialog.getSaveFileName(
            self,
            "Save results as",
            "",
            "Cascade BLAST databases (*.sqlite)"
        )

        if result is None:
            return

        if not result.endswith(".sqlite"):
            result += ".sqlite"

        shutil.copy(self.__db_path, result)
        show_info(
            self,
            "Success",
            f"The results have been save at {result}"
        )

    def __del__(self):
        if self.__tmp_dir is not None:
            self.__tmp_dir.cleanup()

    def __create_initial_table(self):
        organisms = query_organisms(self.__db_path)
        self.__set_organisms(organisms)

    def __toggle_working_widgets(self, visible : bool):

        for widget in self.__action_widgets:
            widget.setEnabled(not visible)

    def __create_cascade(self, args : CreateCascadeDatabaseArgs):

        self.__toggle_working_widgets(True)
        self.__cascade_process = cascade_process = self.__run_create_cascade(args)
        cascade_process.setParent(self)
        cascade_process.finished.connect(self.__on_create_cascade_complete)
        self.__progress_widget.monitor_progress(cascade_process)

    @pyqtSlot(int, QProcess.ExitStatus, name = "__on_create_cascade_complete")
    @with_error_handler()
    def __on_create_cascade_complete(self, code, exit_status = QProcess.NormalExit):

        process = self.__cascade_process
        assert process is not None, "Bug in the code, process is None"

        if process.exitCode() != 0:
            process_output = process.readAllStandardError().data().decode('utf-8')
            message = f"Error running BLAST: {process_output}"
            show_error(self, "Error", message)

        self.__toggle_working_widgets(False)
        self.__create_initial_table()

    def __run_create_cascade(self, args : CreateCascadeDatabaseArgs):

        assert self.__tmp_dir, "Tmp file needed to create a cascade"
        spec_file = path.join(self.__tmp_dir.name, "spec.json")
        with open(spec_file, 'w') as spec_stream:
            args.write_spec(spec_stream)

        fasta_file = path.join(self.__tmp_dir.name, "sequences.fasta")
        SeqIO.write(args.sequences, fasta_file, "fasta")

        return create_cascade(
            self.__db_path,
            spec_file,
            fasta_file,
            args.target_identity,
            args.email,
            args.domain
        )

    def __set_organisms(self, organisms : QueryCascadeResult):
        self.__model = model = OrganismsTableModel(organisms)
        organisms_table = self.__ui.organismsTable
        organisms_table.setModel(model)
        organisms_table.resizeColumnToContents(0)

        self.__policy_combos = {}

        for i,step in enumerate(model.steps):
            filter_widget = CascadeFilterWidget()
            self.__ui.organismsTable.setIndexWidget(
                model.index(0, i + len(model.META_COLS)),
                filter_widget
            )
            self.__policy_combos[step.step_id] = filter_widget
            organisms_table.setColumnWidth(1+i, CascadeFilterWidget.PREFERRED_WIDTH)

        organisms_table.setRowHeight(0, CascadeFilterWidget.PREFERRED_HEIGHT)

    def __filter_results(self):
        filters = dict(
            (step_id, widget.value())
            for step_id, widget in self.__policy_combos.items()
        )

        def accept_organism(organism_entry: QueryCascadeResultOrganismEntry):
            for step in organism_entry.steps:
                identity_filter = filters[step.step_id]

                if identity_filter.operator == CascadeFilterOperator.LessThan \
                    and step.identity >= identity_filter.treshold:
                    return False
                elif identity_filter.operator == CascadeFilterOperator.GreaterThan \
                    and step.identity <= identity_filter.treshold:
                    return False

            return True

        filtered_rows = set(
            row
            for (row, organism_entry) in self.__model.enumerate_organisms_rows()
            if not accept_organism(organism_entry)
        )

        for row,_ in self.__model.enumerate_organisms_rows():
            if row in filtered_rows:
                self.__ui.organismsTable.hideRow(row)
            else:
                self.__ui.organismsTable.showRow(row)
