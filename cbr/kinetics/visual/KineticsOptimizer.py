from PyQt5.QtWidgets import QWidget

from ..data import *
from .Plot import Plot, Point2d, Series
from .Ui_KineticsOptimizer import Ui_KineticsOptimizer

def as_series(runs: KineticsRuns) -> List[Series]:

    return [
        Series(
            name=f"conc {run.run_metadata.concentration}",
            values=[
                Point2d(
                    i * runs.global_attributes.measurement_interval,
                    value
                )
                for i,value in enumerate(run.data)
            ]
        )
        for run in runs.runs
    ]

class KineticsOptimizer(QWidget):

    def __init__(
        self,
        runs: KineticsRuns
    ) -> None:
        super().__init__()

        self.__ui = Ui_KineticsOptimizer()
        self.__ui.setupUi(self)

        self.__plot_widget = Plot()
        self.layout().replaceWidget(
            self.__ui.plotWidget,
            self.__plot_widget
        )

        self.__runs = runs

        self.__render_plots()

    def __render_plots(self):
        exp = as_series(self.__runs)
        self.__plot_widget.set_series(exp)