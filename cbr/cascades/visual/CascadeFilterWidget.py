from enum import Enum
from typing import NamedTuple
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QWidget

from .Ui_CascadeFilterWidget import Ui_CascadeFilterWidget

class CascadeFilterOperator(Enum):
    Any = 0
    GreaterThan = 1
    LessThan = 2


class CascadeFilter(NamedTuple):
    treshold : float
    operator : CascadeFilterOperator

CASCADE_FILTERS = {
    "Any": CascadeFilterOperator.Any,
    "Less": CascadeFilterOperator.LessThan,
    "Greater": CascadeFilterOperator.GreaterThan
}

class CascadeFilterWidget(QWidget):

    PREFERRED_WIDTH = 175
    PREFERRED_HEIGHT = 50

    filter_changed_signal = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.__ui = Ui_CascadeFilterWidget()
        self.__ui.setupUi(self)

        self.__ui.howComboBox.insertItems(0, CASCADE_FILTERS.keys())
        self.__ui.howComboBox.currentTextChanged.connect(self.filter_changed_signal)

        self.__ui.identitySpinBox.valueChanged.connect(self.filter_changed_signal)

    def value(self) -> CascadeFilter:
        treshold = self.__ui.identitySpinBox.value() / 100
        how = CASCADE_FILTERS[self.__ui.howComboBox.currentText()]

        return CascadeFilter(treshold=treshold, operator = how)