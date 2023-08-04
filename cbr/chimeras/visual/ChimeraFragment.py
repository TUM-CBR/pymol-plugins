from PyQt5.QtWidgets import QHeaderView, QWidget
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

    def __on_item_changed(self, item):

        table = self.__ui.fragmentTable
        self.__fragment = Fragment(
            options = list(
                text
                for i in range(self.__ui.fragmentTable.rowCount())
                for item in [table.item(i,0)] if item is not None
                for text in [re.sub("\\s", "", item.text())] if len(text) > 0
            )
        )

    @property
    def fragment(self):
        return self.__fragment