from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QApplication, QHeaderView, QTableWidgetItem, QWidget
import re

from ..Fragment import Fragment

from .Ui_ChimeraFragment import Ui_ChimeraFragment

class ChimeraFragment(QWidget):

    def __init__(self):
        super(ChimeraFragment, self).__init__()
        self.__ui = Ui_ChimeraFragment()
        self.__ui.setupUi(self)

        self.__ui.fragmentTable.adjustSize()
        self.__ui.fragmentTable.itemChanged.connect(self.__on_item_changed)
        
        self.__fragment = Fragment([])

        for _ in range(20):
            self.__ui.fragmentTable.insertRow(0)

        self.__ui.fragmentTable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        self.__ui.pasteButton.clicked.connect(self.paste_from_clipboard)

    def __on_item_changed(self, item):

        table = self.__ui.fragmentTable
        self.__fragment = Fragment(
            options = list(
                text.upper()
                for i in range(self.__ui.fragmentTable.rowCount())
                for item in [table.item(i,0)] if item is not None
                for text in [re.sub("\\s", "", item.text())] if len(text) > 0
            )
        )

    @pyqtSlot()
    def paste_from_clipboard(self):
        tableWidget = self.__ui.fragmentTable
        clipboard = QApplication.clipboard()
        clipboard_text = clipboard.text()

        rows = clipboard_text.split('\n')
        for row_index, row in enumerate(rows):
            columns = row.split('\t')  # Assuming tab-separated values
            tableWidget.insertRow(tableWidget.rowCount())
            for col_index, col_data in enumerate(columns):
                item = QTableWidgetItem(col_data)
                tableWidget.setItem(row_index, col_index, item)

    @property
    def fragment(self):
        return self.__fragment