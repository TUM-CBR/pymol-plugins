from io import StringIO
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTableView, QWidget
from PyQt5.QtGui import QKeyEvent, QGuiApplication, QKeySequence
from typing import Optional, Protocol

class ClipboardHandler(Protocol):

    def on_copy_requested(self, sender: QTableView) -> Optional[str]:
        raise NotImplementedError()

class DefaultClipboardHandler():

    def on_copy_requested(self, sender: QTableView):
        model = sender.model()

        if model is None:
            return None
        
        selected_items = sender.selectedIndexes()

        if len(selected_items) == 0:
            return None
        
        if len(selected_items) == 1:
            selected = selected_items[0]
            return model.data(
                selected,
                Qt.ItemDataRole.DisplayRole
            )
        
        start = min(selected_items, key=lambda i: (i.row(), i.column()))
        end = max(selected_items, key=lambda i: (i.row(), i.column()))

        headers = "\t".join(
            model.headerData(i, Qt.Orientation.Horizontal, Qt.ItemDataRole.DisplayRole)
            for i in range(start.column(), end.column() + 1)
        )

        with StringIO() as buffer:
            #todo ignore items that are not selectred
            buffer.write(headers)
            buffer.write("\n")
            for i in range(start.row(), end.row() + 1):

                buffer.write(
                    "\t".join(
                        (value is not None) and str(value) or " "
                        for j in range(start.column(), end.column() + 1)
                        for index in [model.index(i, j)]
                        for value in [model.data(index, Qt.ItemDataRole.DisplayRole)]
                    )
                )
                buffer.write("\n")
            buffer.seek(0)
            return buffer.read()

class TableViewClipboard(QTableView):

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)

        self.__clipboard_handler : Optional[ClipboardHandler] = DefaultClipboardHandler()

    def setClipboardHandler(self, handler: Optional[ClipboardHandler]) -> None:
        self.__clipboard_handler = handler
    
    def __handle_copy(self) -> bool:

        handler = self.__clipboard_handler
        clipboard = QGuiApplication.clipboard()
        if handler is None or clipboard is None:
            return False
        
        data = handler.on_copy_requested(self)

        if data is None:
            return False

        clipboard.setText(data)

        return True

    def keyPressEvent(self, e: Optional[QKeyEvent]) -> None:

        if e is None:
            return super().keyPressEvent(e)

        if e.matches(QKeySequence.StandardKey.Copy) \
            and self.__handle_copy():
            return

        return super().keyPressEvent(e)