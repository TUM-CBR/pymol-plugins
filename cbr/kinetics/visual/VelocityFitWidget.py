from typing import NamedTuple

from ..compute import ComputeHandler

from .FitWidgetBase import FitWidgetBase
from .Ui_VelocityFitWidget import Ui_VelocityFitWidget

class VelocityConfiguration(NamedTuple):
    samples_to_consider: int

class VelocityFitWidget(FitWidgetBase):

    def __init__(self, compute: ComputeHandler) -> None:
        super().__init__(compute)
        self.__ui = Ui_VelocityFitWidget()

    def __update_busy_state__(self, is_busy: bool):
        self.__ui.fitButton.setEnabled(not is_busy)