from io import StringIO
from typing import Callable
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QAction, QApplication, QMenu, QTableWidget

def open_copy_context_menu(qtable : QTableWidget, pos):

    @pyqtSlot()
    def on_click():
        copy_qtable_to_clipboard(qtable)

    context_menu = QMenu(qtable)
    copy_action = QAction("Copy", qtable)
    copy_action.triggered.connect(on_click)
    context_menu.addAction(copy_action)
    context_menu.exec_(qtable.mapToGlobal(pos))

def get_qtable_content(
    table_widget : QTableWidget,
    is_included : Callable[[int, int], bool]
) -> str:

    if table_widget.rowCount() == 0 \
        or table_widget.columnCount() == 0:
        return ""
    
    with StringIO() as copied_data:
        for row in range(table_widget.rowCount()):
            for column in range(table_widget.columnCount()):
                if is_included(row, column):
                    item = table_widget.item(row, column)
                    assert item is not None, "Item should be found. Bug in the code!"
                    copied_data.write(str(item.text()))
                    copied_data.write("\t")
            copied_data.write("\n")

        copied_data.seek(0)
        return copied_data.read()

def copy_qtable_to_clipboard(table_widget : QTableWidget):
    selected_indexes = table_widget.selectedIndexes()
    if len(selected_indexes) > 0:
        selected_cells = set()
        for index in selected_indexes:
            selected_cells.add((index.row(), index.column()))

        copied_data = get_qtable_content(
            table_widget,
            lambda row, column: (row, column) in selected_cells
        )

        clipboard = QApplication.clipboard()
        clipboard.setText(copied_data)