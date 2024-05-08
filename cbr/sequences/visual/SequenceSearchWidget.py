from io import StringIO
from warnings import warn

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json
from os import path
from PyQt5.QtCore import QModelIndex, QObject, Qt, pyqtSlot, QAbstractTableModel, QThread
from PyQt5.QtWidgets import QDialog, QWidget
import re
from typing import Any, Optional, Sequence

from ...blast.blast import Blast
from ...core import namedtuple
from ...core.Context import Context
from ...core.Qt.QtWidgets import show_error
from ..data import QueryResult, QueryResults, SequencesConfig
from ..process import RunScanResult, RunSearchResult, SequenceCommandRunner
from .OpenOrCreateDialog import OpenOrCreateDialog
from .Ui_SequenceSearchWidget import Ui_SequenceSearchWidget

IS_MAYBE_FASTA_RE = re.compile(r"\s*>.+\n")

class QueryResultModel(QAbstractTableModel):

    COLUMN_ID = "Sequence Id"
    COLUMN_IDENTITY = "Identity"
    COLUMN_FILE = "File"

    COLUMNS = [COLUMN_ID, COLUMN_IDENTITY, COLUMN_FILE]

    def __init__(self, parent: Optional[QObject] = None) -> None:
        super().__init__(parent)

        self.__results: Optional[QueryResults] = None

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        
        if self.__results is None:
            return 0
        else:
            return len(self.__results.results)
        
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.COLUMNS)
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role != Qt.ItemDataRole.DisplayRole:
            return
        
        if orientation == Qt.Orientation.Horizontal:
            return self.COLUMNS[section]

    def __get_record(self, index: QModelIndex) -> QueryResult:

        results = self.__results
        assert results is not None, "Only records that exist should be requested"
        return results.results[index.row()]
    
    def set_result(self, results: QueryResults):
        self.__results = results
        self.modelReset.emit()

    def __display_role(self, index: QModelIndex):

        column_ix = index.column()
        column = self.COLUMNS[column_ix]
        record = self.__get_record(index)

        if column == self.COLUMN_ID:
            return record.id
        elif column == self.COLUMN_IDENTITY:
            return f"{int(100*round(record.identity, 2))}%"
        elif column == self.COLUMN_FILE:
            return record.file_location
        else:
            warn(f"SequenceSearchWidget.py: The column '{column}' is unknown")

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__display_role(index)

class SequenceSearchWidget(QWidget):
    """
    This widget is the main user interface used to search for DNA material that
    matches an input protein sequence.
    """

    def __init__(
        self,
        context: Context,
        parent: Optional[QWidget] = None
    ) -> None:
        super().__init__(parent)
        self.__context = context
        self.__config_directory = context.get_config_directory("sequences")
        self.__config_file = path.join(self.__config_directory, "settings.json")
        self.__ui = Ui_SequenceSearchWidget()
        self.__ui.setupUi(self)
        self.__set_busy(False)

        self.__ui.searchButton.clicked.connect(self.__on_query_clicked)
        self.__ui.rescanButton.clicked.connect(self.__on_rescan_clicked)
        self.__ui.changeDatabase.clicked.connect(self.__on_change_database)

        self.__results_model = QueryResultModel()
        self.__ui.resultsTable.setModel(self.__results_model)

        self.__database: Optional[str] = None
        self.__open_dialog = OpenOrCreateDialog(self.__databases_directory())

        self.__command_thread = QThread(self)
        self.__command_runner = SequenceCommandRunner(Blast(context))
        self.__command_runner.moveToThread(self.__command_thread)
        self.__command_thread.start()
        context.on_app_close(self.__on_app_close)
        self.__command_runner.scan_done_signal.connect(self.__on_scan_done)
        self.__command_runner.query_done_signal.connect(self.__on_query_done)


        # This is always called last in the constructro
        self.__open_default_database()

    def __on_app_close(self):
        self.__command_thread.quit()
        self.__command_thread.wait()

    def __set_busy(self, busy: bool):
        self.__ui.scanProgress.setVisible(busy)
        self.__ui.rescanButton.setEnabled(not busy)
        self.__ui.searchButton.setEnabled(not busy)
        self.__ui.changeDatabase.setEnabled(not busy)

    def __databases_directory(self):
        return self.__config_directory
    
    #def __list_databases(self) -> Sequence[str]:
    #    pattern = path.join(self.__databases_directory(), '*.sqlite')
    #    return list(glob(pattern))
    
    def __write_config(self, config: SequencesConfig):
        with open(self.__config_file, 'w') as config_stream:
            return json.dump(
                config.to_json_dict(),
                config_stream
            )
    
    def __read_config(self) -> SequencesConfig:

        try:
            with open(self.__config_file, 'r') as config:
                return namedtuple.load(SequencesConfig, config)
        except FileNotFoundError:
            return SequencesConfig()
        except ValueError:
            return SequencesConfig()
        
    def __on_scan_done(self, scan_result: RunScanResult):
        self.__set_busy(False)

        if scan_result.error is not None:
            show_error(
                self,
                "Scan Error",
                scan_result.error
            )

    def __scan(self, search_folders: Optional[Sequence[str]] = None):

        database = self.__database

        if database is None:
            show_error(self, "Scan Error", "No database has been selected!")
            return

        self.__set_busy(True)
        self.__command_runner.run_scan(database, search_folders)

    @pyqtSlot()
    def __on_rescan_clicked(self):
        self.__scan()

    def __get_seqs(self) -> Optional[Sequence[SeqRecord]]:

        seqs_text = self.__ui.sequenceText.toPlainText()

        if IS_MAYBE_FASTA_RE.search(seqs_text) is None:
            seqs_text = f">query\n{seqs_text}"

        try:
            return list(SeqIO.parse(StringIO(seqs_text), format='fasta'))
        except ValueError:
            return None
        
    def __on_query_done(self, result: RunSearchResult):

        self.__set_busy(False)

        if result.error is not None:
            show_error(
                self,
                "Query Error",
                "Failed to query the database: " + result.error
            )
        else:
            assert result.results is not None, "Results should not be None if there is no error."
            self.__results_model.set_result(result.results)
            self.__ui.resultsTable.resizeColumnsToContents()

    def __query(self):

        database = self.__database

        if database is None:
            show_error(self, "Query Error", "No database has been selected")
            return
        
        seqs = self.__get_seqs()

        if seqs is None or len(seqs) < 1:
            show_error(
                self,
                "Query Error",
                "The given query is not a valid sequence."
            )
            return

        self.__set_busy(True)
        self.__command_runner.run_query(database, seqs)

    @pyqtSlot()
    def __on_query_clicked(self):
        self.__query()

    def __open_or_create_database(self) -> bool:

        if self.__open_dialog.exec() == QDialog.DialogCode.Accepted:
            value = self.__open_dialog.value()
            assert value is not None, "Bug in OpenOrCreateDialog, returned None upon success."
            self.__database = value.db_file
            self.__scan(value.search_folders)
            config = self.__read_config()
            self.__write_config(
                config._replace(selected_database = value.db_file)
            )
            return True
        else:
            return False
        
    @pyqtSlot()
    def __on_change_database(self):
        self.__open_or_create_database()

    def __open_default_database(self):

        config = self.__read_config()

        if config.selected_database is not None:
            self.__database = config.selected_database
        elif not self.__open_or_create_database():
            self.__context.close_app()