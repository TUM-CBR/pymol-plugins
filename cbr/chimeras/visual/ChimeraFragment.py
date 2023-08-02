from PyQt5.QtWidgets import QWidget

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

    def __on_item_changed(self, item):

        table = self.__ui.fragmentTable
        self.__fragment = Fragment(
            options = list(
                 item.text()
                for i in range(self.__ui.fragmentTable.rowCount())
                for item in [table.itemAt(i, 0)] if item is not None
            )
        )

    @property
    def fragment(self):
        return self.__fragment