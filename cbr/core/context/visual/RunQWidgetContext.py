from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QCloseEvent
from PyQt5.QtWidgets import QMainWindow, QWidget
from typing import Callable, Optional

from .Ui_RunQWidgetContext import Ui_RunQWidgetContext

class RunQWidgetContext(QMainWindow):
    def __init__(
        self,
        widget : QWidget,
        on_close : 'Callable[[RunQWidgetContext], None]',
        user_manual: QUrl
    ):
        super().__init__()
        self.__ui = Ui_RunQWidgetContext()
        self.__ui.setupUi(self)
        self.__ui.userManualLabel.setText(
            MANUAL_TEXT_TEMPLATE.format(
                user_manual=user_manual.url()
            )
        )

        self.__on_close = on_close
        self.setAppWidget(widget)

    def setAppWidget(self, widget: QWidget) -> None:

        layout = self.__ui.mainWindow.layout()
        assert layout is not None, "There must be a layout"
        layout.replaceWidget(
            self.__ui.centralWidget,
            widget
        )
        self.__ui.centralWidget = widget

    def closeEvent(self, a0: Optional[QCloseEvent]) -> None:
        super().closeEvent(a0)
        self.__on_close(self)

MANUAL_TEXT_TEMPLATE = \
"""
    <html>
        <head/>
            <body>
            <p><a href="{user_manual}"><span style=" text-decoration: underline; color:#0000ff;">User Manual</span></a></p>
            </body>
        </html>
"""