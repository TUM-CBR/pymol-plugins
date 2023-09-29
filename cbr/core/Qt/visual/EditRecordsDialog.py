from typing import Any, Callable, Dict, Generic, NamedTuple, Optional, Type, TypeVar
from PyQt5.QtCore import QModelIndex, QObject, Qt, pyqtSlot
from PyQt5.QtGui import QDoubleValidator, QIntValidator, QValidator
from PyQt5.QtWidgets import QDialog, QItemDelegate, QLineEdit, QStyleOptionViewItem, QTableWidgetItem,QWidget

from .Ui_editrecords import Ui_EditRecords

K_ARGS_META = "__args_meta__"
VALUE_CONTEXT_ROLE = Qt.UserRole

COLUMN_NAMES = ["Parameter", "Value", "Description"]
NUM_COLUMNS = len(COLUMN_NAMES)

T = TypeVar('T', bound=NamedTuple)
CellDataParser = Callable[[str], Any]

class SupportedType(NamedTuple):
    type_class : Type
    parser : CellDataParser
    validator : QValidator
    default : Any

    def is_supported(self, other_type : Type) -> bool:
        return self.type_class == other_type

SUPPORTED_NATIVE_TYPES = [
    SupportedType(type_class=int, parser=lambda x: int(x), validator=QIntValidator(), default=0),
    SupportedType(type_class=float, parser=lambda x: float(x), validator=QDoubleValidator(), default=0)
]

def argsinfo(**args):

    def decorator(structure : Type[T]) -> Type[T]:

        if not issubclass(structure, tuple):
            raise ValueError("The decorator 'argsinfo' can only with used with subclasses of 'NamedTuple'")

        for key in args.keys():
            if key not in structure._fields:
                raise ValueError(f"Argument {key} does not exist in {structure.__class__}")

        setattr(structure, K_ARGS_META, args)

        return structure

    return decorator

class ValueContext(NamedTuple):
    field_name : str
    field_parser : Optional[Callable[[str], Any]]

    def parse(self, text : str) -> Any:
        if self.field_parser is None:
            return text

        return self.field_parser(text)

class DelegateWithValidator(QItemDelegate):

    def __init__(self, parent: QObject, validator : QValidator):
        super().__init__(parent)
        self.__validator = validator

    def createEditor(
        self,
        parent: QWidget,
        option: QStyleOptionViewItem,
        index: QModelIndex
    ) -> QWidget:

        editor = QLineEdit(parent)
        editor.setValidator(self.__validator)
        return editor


class EditRecordsDialog(QDialog, Generic[T]):

    def __init__(
        self,
        parent: QWidget,
        current : T,
        defaults : T
    ) -> None:
        super(EditRecordsDialog, self).__init__(parent)

        self.__ui = Ui_EditRecords()
        self.__ui.setupUi(self)
        self.__current = current
        self.__defaults = defaults
        self.__ui.valuesTable.cellChanged.connect(self.__on_item_changed)
        self.__ui.saveButton.clicked.connect(self.__on_close)
        self.__ui.restoreDefaultsButton.clicked.connect(self.__on_restore_defaults)

        #render the table
        self.__render()

    def __validate(self):
        
        if list(self.__current._field_types.items()) != list(self.__defaults._field_types.items()):
            raise ValueError("Default value must have the same fields as current value.")

    @pyqtSlot()
    def __on_restore_defaults(self):
        self.__current = self.__defaults
        self.__render()

    @pyqtSlot()
    def __on_close(self):
        self.accept()

    @pyqtSlot(int, int)
    def __on_item_changed(self, row: int, col: int):

        if col != 1:
            return

        item = self.__ui.valuesTable.item(row, col)
        if item is None:
            return

        field_context : ValueContext = item.data(VALUE_CONTEXT_ROLE)
        field = field_context.field_name
        value = field_context.parse(item.text())

        self.__current = self.__current._replace(**{field : value})

    def new_value(self) -> T:
        return self.__current

    def __render(self):

        self.__validate()

        current = self.__current
        fields = list(current._field_types.items())
        meta : Dict[str, str] = getattr(current, K_ARGS_META, {})

        table = self.__ui.valuesTable
        table.clear()
        table.setRowCount(len(fields))
        table.setColumnCount(NUM_COLUMNS)
        table.setHorizontalHeaderLabels(COLUMN_NAMES)

        row = 0
        for (name, value_type) in fields:
            description = meta.get(name)
            value = None
            delegate = None
            parser : Optional[CellDataParser] = None

            for supported_type in SUPPORTED_NATIVE_TYPES:
                if supported_type.is_supported(float):
                    value = getattr(current, name, supported_type.default)
                    delegate = table.setItemDelegateForRow(
                        row,
                        DelegateWithValidator(table, supported_type.validator)
                    )
                    parser = supported_type.parser

            if value is not None:

                # The cell displaying the name of the attribute
                name_cell = QTableWidgetItem(name)
                name_cell.setFlags(name_cell.flags() & ~Qt.ItemIsEditable)
                table.setItem(row, 0, name_cell)

                # The cell displaying the value of the attribute
                value_cell = QTableWidgetItem(str(value))
                value_cell.setData(VALUE_CONTEXT_ROLE, ValueContext(field_name=name, field_parser=parser))
                table.setItem(row, 1, value_cell)
                if delegate is not None:
                    table.setItemDelegateForRow(row, delegate)

                # The cell displaying the description of the attribute
                description_cell = QTableWidgetItem(description or "")
                description_cell.setFlags(name_cell.flags() & ~Qt.ItemIsEditable)
                table.setItem(row, 2, description_cell)

                row += 1

        table.resizeColumnsToContents()