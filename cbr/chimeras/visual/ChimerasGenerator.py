import itertools
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QHeaderView, QListWidgetItem, QTableWidgetItem, QWidget
from typing import Dict

from ...core.Context import Context

from .ChimeraFragment import ChimeraFragment
from .Ui_ChimerasGenerator import Ui_ChimerasGenerator

class ChimerasGenerator(QWidget):

    def __init__(self, context : Context):

        super(ChimerasGenerator, self).__init__()
        self.__ui = Ui_ChimerasGenerator()
        self.__ui.setupUi(self)
        self.__editors : Dict[int, ChimeraFragment] = {}
        self.__set_fragments_count(10)

        self.__ui.fragmentsList.itemClicked.connect(self.__selected_item)
        self.__ui.fragmentsStack.setCurrentIndex(0)
        self.__ui.generateButton.clicked.connect(self.__generate_sequences)
        self.__ui.resultsTable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

    @pyqtSlot(QListWidgetItem)
    def __selected_item(self, item : QListWidgetItem):
        index = self.__ui.fragmentsList.indexFromItem(item)
        self.__ui.fragmentsStack.setCurrentIndex(index.row())

    def __set_fragments_count(self, count : int):

        fragmentsList = self.__ui.fragmentsList
        fragmentsStack = self.__ui.fragmentsStack
        current_count = fragmentsList.count()
        items_to_add = max(count - current_count, 0)
        items_to_rm = max(current_count - count, 0)

        for i in range(current_count, current_count + items_to_add):
            number = i + 1
            fragmentsList.insertItem(i, "fragment %i" % number)
            self.__editors[i] = ChimeraFragment()
            fragmentsStack.insertWidget(i, self.__editors[i])

        for number in range(items_to_rm):
            current_count = fragmentsList.count() - 1
            fragmentsList.takeItem(current_count)
            fragmentsStack.removeWidget(self.__editors[current_count])
            del self.__editors[current_count]

    def __generate_sequences(self):
        catalogue = (
            fragment.fragment.options
            for fragment in self.__editors.values()
            if len(fragment.fragment.options) > 0
        )

        resutlsTable = self.__ui.resultsTable
        resutlsTable.clearContents()
        results = list(enumerate(itertools.product(*catalogue)))
        resutlsTable.setRowCount(len(results))

        for i,item in results:
            row = QTableWidgetItem("".join(item))
            resutlsTable.setItem(i, 0, row)
