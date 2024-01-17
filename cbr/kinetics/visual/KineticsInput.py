from PyQt5.QtCore import pyqtSignal, pyqtSlot, Qt
from PyQt5.QtWidgets import QApplication, QTableWidgetItem, QWidget
from typing import Iterable, cast

from ...core.Qt.visual.NamedTupleEditor import namedtuple_eidtor

from ..data import *
from .Ui_KineticsInput import Ui_KineticsInput

NUM_SERIES = 50

class KineticsInput(QWidget):

    runs_selected = pyqtSignal()

    def __init__(self) -> None:
        super().__init__()
        self.__ui = Ui_KineticsInput()
        self.__ui.setupUi(self)

        self.__globals_editor = namedtuple_eidtor(
            self.__ui.runGlobalsTable,
            GlobalAttributes(3125, 1, 10, 0.001)
        )

        self.__runs_metadata_editor = namedtuple_eidtor(
            self.__ui.runMetadataTable,
            *[None for _ in range(0,NUM_SERIES)],
            tuple_type=RunMetadata,
            default_item=RunMetadata()
        )

        self.__ui.runDataTable.setRowCount(100)
        self.__ui.runDataTable.setColumnCount(NUM_SERIES)

        self.__ui.pasteParametersButton.clicked.connect(self.__on_paste_parameters)
        self.__ui.pasteDataButton.clicked.connect(self.__on_paste_data)
        self.__ui.proceedButton.clicked.connect(self.runs_selected)

    @pyqtSlot()
    def __on_paste_parameters(self):
        clipboard = QApplication.clipboard()

        assert clipboard is not None, "This operation requires a GUI"

        self.__runs_metadata_editor.write_string(clipboard.text())

    @pyqtSlot()
    def __on_paste_data(self):
        clipboard = QApplication.clipboard()

        assert clipboard is not None, "This operation requires a GUI"

        for row, line in enumerate(clipboard.text().splitlines()):
            for col, data in enumerate(line.split('\t')):
                self.__ui.runDataTable.setItem(
                    row,
                    col,
                    QTableWidgetItem(data)
                )

    def __get_series(self) -> Iterable[KineticsRun]:
        runs = self.__runs_metadata_editor.current_values

        for (i, run) in enumerate(runs):

            if run is None:
                break

            data = [
                float(value.data(Qt.ItemDataRole.DisplayRole))
                for j in range(0, self.__ui.runDataTable.rowCount())
                for value in [self.__ui.runDataTable.item(j, i)] if value is not None
            ]

            yield KineticsRun(
                run_metadata = run,
                data = data
            )

    def get_runs(self) -> KineticsRuns:
        global_attributes = self.__globals_editor.current_values[0]

        return KineticsRuns(
            global_attributes = cast(GlobalAttributes, global_attributes),
            runs = list(self.__get_series())
        )