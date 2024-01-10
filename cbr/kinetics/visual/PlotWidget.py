from PyQt5.QtWidgets import QWidget

from .Plot import Plot, SeriesSet
from .SeriesModel import SeriesModel
from .Ui_PlotWidget import Ui_PlotWidget

class PlotWidget(QWidget):

    def __init__(self) -> None:
        super().__init__()

        self.__ui = Ui_PlotWidget()
        self.__ui.setupUi(self)

        self.__plot = Plot()

        plot_layout = self.__ui.plotContainer.layout()

        assert plot_layout is not None, "Bug in the code. No layout for plot container."

        plot_layout.replaceWidget(
            self.__ui.plotWidget,
            self.__plot
        )

    def set_series(self, series: SeriesSet):
        self.__plot.set_series(series)
        self.__ui.dataTable.setModel(
            SeriesModel(
                series.series,
                render_series_header=lambda i,s: s.metadata.name
            )
        )

    