from PyQt5.QtCore import QStringListModel
from PyQt5.QtWidgets import QWidget

from ...core.Context import Context

from .Ui_ChimerasGenerator import Ui_ChimerasGenerator

class ChimerasGenerator(QWidget):

    def __init__(self, context : Context):

        super(ChimerasGenerator, self).__init__()
        self.__ui = Ui_ChimerasGenerator()
        self.__ui.setupUi(self)
        self.__fragments_count = 0
        self.__set_fragments_count(5)

    def __set_fragments_count(self, count : int):

        fragmentsList = self.__ui.fragmentsList
        current_count = fragmentsList.count()
        items_to_add = max(count - current_count, 0)
        items_to_rm = max(current_count - count, 0)

        for i in range(current_count, current_count + items_to_add):
            fragmentsList.insertItem(i, "fragment %i" % i)



