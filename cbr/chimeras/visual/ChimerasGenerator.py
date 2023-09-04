from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import QIntValidator
from PyQt5.QtWidgets import QHeaderView, QListWidgetItem, QMessageBox, QWidget
import re
from typing import Dict, Iterable

from ...core.Context import Context

from .ChimeraFragment import ChimeraFragment
from .ChimerasResults import ChimerasGeneratorArgs, ChimeraGeneratorPosition, ChimerasResults
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
        self.__results_widget = None

        self.__ui.overhangsTable.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.__ui.overhangsTable.setRowCount(20)
        self.__ui.overhangLenght.setValidator(QIntValidator(0,10))

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
            fragmentsList.insertItem(i, "position %i" % number)
            self.__editors[i] = ChimeraFragment()
            fragmentsStack.insertWidget(i, self.__editors[i])

        for number in range(items_to_rm):
            current_count = fragmentsList.count() - 1
            fragmentsList.takeItem(current_count)
            fragmentsStack.removeWidget(self.__editors[current_count])
            del self.__editors[current_count]

    def __get_binding_codons(self) -> Iterable[str]:

        overhangs_table = self.__ui.overhangsTable
        return(
            overhang
            for i in range(overhangs_table.rowCount())
            for item in [overhangs_table.item(i, 0)] if item is not None
            for overhang in [re.sub("\\s", "", item.text())] if len(overhang) > 0
        )

    def __generate_sequences(self):

        positions = [
            ChimeraGeneratorPosition(sequences = fragment.fragment.options)
            for fragment in self.__editors.values()
            if len(fragment.fragment.options) > 0
        ]

        try:
            self.__results_widget = widget = ChimerasResults(
                ChimerasGeneratorArgs.generate_chimeras(
                    positions,
                    self.__get_binding_codons(),
                    int(self.__ui.overhangLenght.text())
                )
            )
            widget.show()
        except Exception as e:
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Critical)
            msgBox.setText(str(e))
            msgBox.setWindowTitle("Error")
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec()