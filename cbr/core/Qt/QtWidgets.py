from concurrent.futures import Future
import csv
from io import StringIO
from typing import Any, Callable, Iterable, Optional, Set, TextIO, TypeVar, cast
from PyQt5.QtCore import QTimer, pyqtSignal, pyqtSlot, QObject, Qt
from PyQt5.QtWidgets import QAction, QApplication, QLabel, QMenu, QMessageBox, QProgressBar, QSlider, QTableWidget, QTableView, QWidget

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
    is_included : Callable[[int, int], bool],
    column_headers : bool = True
) -> str:

    if table_widget.rowCount() == 0 \
        or table_widget.columnCount() == 0:
        return ""

    columns : Set[int] = set()
    
    with StringIO() as copied_data:
        for row in range(table_widget.rowCount()):
            for column in range(table_widget.columnCount()):
                if is_included(row, column):

                    if column_headers:
                        columns.add(column)
                    item = table_widget.item(row, column)
                    assert item is not None, "Item should be found. Bug in the code!"
                    copied_data.write(str(item.text()))
                    copied_data.write("\t")
            copied_data.write("\n")

        columns_text = "\t".join(
            table_widget.horizontalHeaderItem(i).text() # type: ignore[reportOptionalMemberAccess]
            for i in sorted(columns)
        )

        copied_data.seek(0)
        return columns_text + "\n" + copied_data.read()

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

def show_info(
    parent : Optional[QWidget],
    title : str,
    description : str = '',
    window_title : str = 'Information'
):
    dialog = QMessageBox(parent)
    dialog.setIcon(QMessageBox.Information)
    dialog.setWindowTitle(window_title)
    dialog.setText(title)
    dialog.setInformativeText(description)
    dialog.setStandardButtons(QMessageBox.Ok)
    dialog.exec_()

def show_error(
    parent : Optional[QWidget],
    title : str,
    description : str = '',
    window_title : str = 'Error'
):
    error_dialog = QMessageBox(parent)
    error_dialog.setIcon(QMessageBox.Critical)
    error_dialog.setWindowTitle(window_title)
    error_dialog.setText(title)
    error_dialog.setInformativeText(description)
    error_dialog.setStandardButtons(QMessageBox.Ok)
    error_dialog.exec_()

def show_exception(
    parent : QWidget,
    exn : BaseException
):
    return show_error(
        parent,
        title=exn.__class__.__name__,
        description=str(exn)
    )

TResult = TypeVar('TResult')

TErrorHandler = TypeVar('TErrorHandler')

class ThrottleHandler(QObject):

    def __init__(
        self,
        time : int,
        widget : QObject,
        action
    ):
        super().__init__()
        self.__action = action
        self.__last_args = []
        self.__last_kwargs = {}
        self.__widget = widget
        timer = self.__timer = QTimer()
        timer.setInterval(time)
        timer.setSingleShot(True)
        timer.timeout.connect(self.__on_timeout)

    def __init_timer__needed(self):

        if self.__timer.isActive():
            return

        self.__timer.start()

    def __call__(self, *args, **kwargs):

        self.__init_timer__needed()
        self.__last_args = args
        self.__last_kwargs = kwargs

    @pyqtSlot()
    def __on_timeout(self):
        self.__action(self.__widget, *self.__last_args, **self.__last_kwargs)


class ThrottleFactory:

    def __init__(self, time, fn):
        self.__time = time
        self.__fn = fn
        self.__handlers = {}

    @pyqtSlot(QObject)
    def __on_widget_destroyed(self, widget: QWidget):
        self.__handlers.pop(widget)

    def __call__(self, widget, *args, **kwargs):

        handler = self.__handlers.get(widget)

        if handler is None:
            handler = self.__handlers[widget] = ThrottleHandler(self.__time, widget, self.__fn)
            #widget.destroyed.connect(self.__on_widget_destroyed)

        return handler(*args, **kwargs)

def throttle(time : int):

    def handler(fn : TErrorHandler) -> TErrorHandler:

        factory = ThrottleFactory(time, fn)

        def run(widget, *args, **kwargs):
            factory(widget, *args, **kwargs)
        return cast(TErrorHandler, run)

    return handler

def with_error_handler(*args, **kwargs):

    invalid_use_exception = ValueError("'with_error_handler' can only be used with methods of QWidget instances!")

    def function_factory(fn : TErrorHandler) -> TErrorHandler:

        def new_function(qWidget : QWidget, *args, **oargs):

            nonlocal kwargs
            if not isinstance(qWidget, QWidget):
                raise invalid_use_exception

            result = None
            try:
                result = cast(Any, fn)(qWidget, *args, **oargs)
            except Exception as exn:
                show_exception(qWidget, exn)
                return

            if result is not None:
                raise ValueError("'with_error_handler' can only be used for methods that return None")

        return cast(TErrorHandler, new_function)

    return function_factory

class ProgressManager(QObject):

    on_result = pyqtSignal(object)
    on_exception = pyqtSignal(Exception)

    def __init__(
        self,
        progress : QProgressBar,
        disable : Iterable[QWidget] = []
    ):

        super(ProgressManager, self).__init__()
        self.__progress = progress
        self.__disable = list(disable)
        self.__set_ready()

    def __set_ready(self):
        self.__progress.setVisible(False)

        for widget in self.__disable:
            widget.setEnabled(True)

    def __set_busy(self):
        self.__progress.reset()
        self.__progress.setMinimum(0)
        self.__progress.setMaximum(0)
        self.__progress.setVisible(True)

        for widget in self.__disable:
            widget.setEnabled(False)

    def __notify_future(self, future : 'Future[Any]') -> bool:

        if not future.done():
            return False

        exception = future.exception()
        if exception:
            self.on_exception.emit(exception)
        else:
            self.on_result.emit(future.result())

        return True

    def watch_progress(self, result : 'Future[Any]'):

        if self.__notify_future(result):
            return

        def __done__(_: Any):
            self.__set_ready()
            self.__notify_future(result)

        result.add_done_callback(__done__)
        self.__set_busy()

    def with_default_error_handler(self, widget : QWidget):

        @pyqtSlot()
        def __handle_error__(exn : Exception):
            show_exception(widget, exn)

        self.on_exception.connect(__handle_error__)

def progress_manager(
    progress : QProgressBar,
    *disable : QWidget
):
    return ProgressManager(progress, disable)

def export_table_view_to_csv(table_view : QTableView, text_stream: TextIO) -> None:
    """
    Export data from a QTableView to a CSV through a text stream.

    :param table_view: QTableView object containing data to be exported.
    :param text_stream: TextIO, a text stream to write the CSV data.
    """
    model = table_view.model()  # QAbstractItemModel
    if not model:
        raise ValueError("The table view does not have a model.")

    writer = csv.writer(text_stream)
    
    # Write header
    headers = []
    for col in range(model.columnCount()):
        headers.append(model.headerData(col, orientation=Qt.Horizontal))
    writer.writerow(headers)
    
    # Write data
    for row in range(model.rowCount()):
        row_data = []
        for col in range(model.columnCount()):
            index = model.index(row, col)
            row_data.append(model.data(index))
        writer.writerow(row_data)

class SliderWithLabel(QObject):

    value_changed = pyqtSignal()
    slider_released = pyqtSignal()

    def __init__(
        self,
        slider : QSlider,
        label : QLabel,
        factor : float = 1
    ):
        super().__init__()

        self.__label = label
        self.__slider = slider
        self.__factor = factor
        label.setText(self.__format(self.value))

        self.__slider.valueChanged.connect(self.__on_value_changed)
        self.__slider.sliderReleased.connect(self.slider_released)

    def reset_label(self):
        value = self.value
        self.__label.setText(self.__format(value))

    @property
    def value(self):
        return self.__slider.value() * self.__factor

    def __format(self, value : float):
        return str(round(value, 2))

    @pyqtSlot()
    def __on_value_changed(self):
        self.reset_label()
        self.value_changed.emit()

slider_with_label = SliderWithLabel