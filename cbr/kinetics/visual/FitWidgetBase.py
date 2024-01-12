from PyQt5.QtCore import Qt, pyqtSlot
from PyQt5.QtWidgets import QWidget

from ..compute import ComputeHandler
from ..data import KineticsRuns

class FitWidgetBase(QWidget):

    def __init__(self, compute: ComputeHandler) -> None:
        super().__init__()
        self.__compute = compute
        self.__compute.on_busy_changed.connect(self.__on_busy_changed)

    def __compute__(self):
        return self.__compute
    
    @pyqtSlot()
    def __on_busy_changed(self):
        is_busy = self.__compute.is_busy()
        self.__update_busy_state__(is_busy)
    
    def __update_busy_state__(self, is_busy: bool):
        pass

    def on_runs_updated(self, runs: KineticsRuns):
        pass