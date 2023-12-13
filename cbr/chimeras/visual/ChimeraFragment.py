from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QHeaderView, QTableWidgetItem, QWidget
import re

from ..Fragment import Fragment

from .Ui_ChimeraFragment import Ui_ChimeraFragment

class ChimeraFragment(QWidget):

    IX_NAME = 0
    IX_SEQ = 1
    ROWS = 20

    def __init__(self):
        super(ChimeraFragment, self).__init__()
        self.__ui = Ui_ChimeraFragment()
        self.__ui.setupUi(self)

        self.__ui.fragmentTable.adjustSize()
        self.__ui.fragmentTable.itemChanged.connect(self.__on_item_changed)
        
        self.__fragment = Fragment(options={})

        for row in range(20):
            self.__ui.fragmentTable.insertRow(0)
            self.__ui.fragmentTable.setItem(0, self.IX_NAME, QTableWidgetItem(str(20 - row)))

        header = self.__ui.fragmentTable.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        header.setSectionResizeMode(self.IX_NAME, QHeaderView.Fixed)
        header.resizeSection(self.IX_NAME, 50)

        self.__ui.pasteButton.clicked.connect(self.paste_from_clipboard)

    def __on_item_changed(self, item):

        table = self.__ui.fragmentTable
        self.__fragment = Fragment(
            options = dict(
                (name, text.upper())
                for i in range(self.__ui.fragmentTable.rowCount())
                for (i_name, i_text) in [(table.item(i, self.IX_NAME), table.item(i, self.IX_SEQ))] \
                    if i_name is not None and i_text is not None
                for (name, text) in [(i_name.text(), re.sub("\\s", "", i_text.text()))] if len(text) > 0
            )
        )

    @pyqtSlot()
    def paste_from_clipboard(self):
        tableWidget = self.__ui.fragmentTable
        clipboard = QApplication.clipboard()
        clipboard_text = clipboard.text()

        selected = tableWidget.selectedIndexes()

        if len(selected) == 1:
            item = selected[0]
            row_offset = item.row()
            col_offset = item.column()
        else:
            row_offset = 0
            col_offset = 0

        row_count = tableWidget.rowCount() - row_offset
        col_count = tableWidget.columnCount() - col_offset

        rows = clipboard_text.split('\n')[:row_count]
        for row_index, row in enumerate(rows):
            columns = row.split('\t')[:col_count]  # Assuming tab-separated values
            for col_index, col_data in enumerate(columns):
                if col_data.strip() == "":
                    continue
                item = QTableWidgetItem(col_data)
                tableWidget.setItem(row_index + row_offset, col_index + col_offset, item)

    @property
    def fragment(self):
        return self.__fragment