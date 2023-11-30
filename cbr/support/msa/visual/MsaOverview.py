from typing import Optional, Tuple
from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import QPoint, Qt, pyqtSignal
from PyQt5.QtGui import QColor, QMouseEvent, QPainter, QPaintEvent, QPixmap, QTransform
from PyQt5.QtWidgets import QWidget

from ....core.Qt.QtGui import paint
from ....clustal import msa
from .MsaMask import MsaMask, MaskPositionMode, RESIDUE_COLORS
from .Ui_MsaOverview import Ui_MsaOverview

MsaSelectedRange = Tuple[int, int]

class MsaOverviewCanvas(QWidget):

    range_selected = pyqtSignal()

    def __init__(self, alignment : MultipleSeqAlignment):
        super().__init__()
        self.__alignment = alignment
        self.__mask = MsaMask.default_mask(alignment)
        self.__image = self.__make_pixmap()
        self.setMouseTracking(False)
        self.__selection = None
        self.__selected_range = None
        self.__render_transform = QTransform().scale(5, 1)

    def mouseMoveEvent(self, a0: QMouseEvent) -> None:

        prev = self.__selection or (a0.pos(), a0.pos())

        self.__selection = (prev[0], a0.pos())
        self.repaint()

    def mousePressEvent(self, a0: QMouseEvent) -> None:
        self.__selection = (a0.pos(), a0.pos())

    @property
    def selected_range(self) -> Optional[MsaSelectedRange]:
        return self.__selected_range

    def mouseReleaseEvent(self, a0: QMouseEvent) -> None:

        if self.__selection is None:
            return

        start, end = self.__selection
        inverted, success = self.__render_transform.inverted()
        assert(success)
        start, end = inverted.map(start), inverted.map(end)

        self.__selection = None

        start_x = start.x()
        end_x = end.x()

        _, msa_start = self.__get_position(0, start_x)
        _, msa_end = self.__get_position(0, end_x)
        self.__selected_range = (msa_start, msa_end)

        self.range_selected.emit()
        self.repaint()

    def set_mask(self, mask: MsaMask):
        self.__mask = mask
        self.__image = self.__make_pixmap()
        self.repaint()

    def __get_alignment_size(self):
        mask = self.__mask
        if mask.mask_mode == MaskPositionMode.HIDE:
            return (mask.row_count, mask.col_count)
        else:
            return (len(self.__alignment), self.__alignment.get_alignment_length())

    def __get_position(self, row: int, col: int) -> Tuple[int, int]:
        mask = self.__mask
        if mask.mask_mode == MaskPositionMode.HIDE:
            return mask.get_position(row, col)
        else:
            return (row, col)

    def __make_pixmap(self) -> QPixmap:

        mask = self.__mask
        rows, cols = self.__get_alignment_size()
        pixmap = QPixmap(cols, rows)

        with paint(pixmap) as manager:
            painter = manager.painter

            for i in range(0, rows):
                for j in range(0, cols):
                    row, col = self.__get_position(i, j)
                    resi = self.__alignment[row][col]
                    assert isinstance(resi, str)
                    if mask.is_masked(row, col):
                        painter.setPen(Qt.black)
                    elif msa.is_blank(resi):
                        painter.setPen(Qt.white)
                    else:
                        painter.setPen(RESIDUE_COLORS[resi])

                    painter.drawPoint(j, i)

        return pixmap

    def __paint_selection(self, painter: QPainter):

        selection = self.__selection
        if selection is None:
            return

        startPoint, endPoint = selection
        startX = startPoint.x()
        endX = endPoint.x()

        if startX > endX:
            tmp = startX
            startX = endX
            endX = tmp

        painter.fillRect(
            startX,
            0,
            endX - startX,
            self.__image.height(),
            QColor(0, 0, 255, 128)
        )

    def paintEvent(self, a0: QPaintEvent) -> None:

        super().paintEvent(a0)

        with paint(self) as manager:
            pximap_size = self.__image.size()
            size = self.__render_transform.map(
                QPoint(pximap_size.width(), pximap_size.height())
            )
            painter = manager.painter
            painter.setTransform(self.__render_transform)
            self.setFixedSize(size.x(), size.y())
            painter.drawPixmap(0,0,self.__image)

        with paint(self) as manager:
            self.__paint_selection(manager.painter)    

class MsaOverview(QWidget):

    range_selected = pyqtSignal()

    def __init__(self, alignment : MultipleSeqAlignment):
        super().__init__()

        self.__ui = Ui_MsaOverview()
        self.__ui.setupUi(self)
        self.__msa_canvas = MsaOverviewCanvas(alignment)
        self.__ui.msaOverviewScrollArea.setWidget(self.__msa_canvas)
        self.__msa_canvas.range_selected.connect(self.range_selected)

    def set_msa_mask(self, mask : MsaMask):
        self.__msa_canvas.set_mask(mask)

    @property
    def selected_range(self) -> Optional[MsaSelectedRange]:
        return self.__msa_canvas.selected_range