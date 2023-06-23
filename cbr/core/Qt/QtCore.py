from PyQt5.QtCore import QObject, pyqtSlot, QTimer

class Throttle(QObject):

    def __init__(
        self,
        timeout : int,
        action,
        *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__timer = QTimer()
        self.__timer.setSingleShot(True)
        self.__timer.timeout.connect(self.__on_timeout)
        self.__timeout = timeout
        self.__action = action
        self.__args = []
        self.__kwargs = {}

    @pyqtSlot()
    def __on_timeout(self):
        self.__action(*self.__args, **self.__kwargs)

    def trigger(self, *args, **kwargs):
        self.__args = args
        self.__kwargs = kwargs
        self.__timer.start(self.__timeout)