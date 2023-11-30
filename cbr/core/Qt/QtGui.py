from PyQt5.QtGui import QPainter, QPaintDevice

class QPainterContextManager:

    def __init__(self, widget: QPaintDevice):
        self.__painter = QPainter()
        self.__widget = widget

    @property
    def widget(self):
        return self.__widget

    @property
    def painter(self):
        return self.__painter

    def __enter__(self):
        self.__painter.begin(self.__widget)

        return self

    def __exit__(self, *args, **kwargs):
        self.__painter.end()

def paint(widget : QPaintDevice) -> QPainterContextManager:
    return QPainterContextManager(widget)