from PyQt5.QtCore import QObject
from PyQt5.QtWidgets import QComboBox


from ..blosum import blosum62, blosum80

SUBSTITUTION_OPTIONS = {
    "SHCEMA classic": None,
    "BLOSUM 80": blosum80,
    "BLOSUM 62": blosum62
}

class SubstitutionSelector(QObject):

    def __init__(self, comboBox : QComboBox):
        super(SubstitutionSelector, self).__init__()
        self.__combo_box = comboBox

        self.__combo_box.addItems(SUBSTITUTION_OPTIONS)

    @property
    def selection(self):
        return SUBSTITUTION_OPTIONS[self.__combo_box.currentText()]