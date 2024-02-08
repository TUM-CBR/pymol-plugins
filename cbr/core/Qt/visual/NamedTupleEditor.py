from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtWidgets import QTableView
from typing import Any, Callable, Iterable, cast, Dict, List, Generic, NamedTuple, Optional, Type, TypeVar

from numpy import append

from ..QtWidgets import show_exception
from .value_handlers import *

OVERRIDES_FIELD = '__editor_field_overrides__'

class MetaFieldOverrides(NamedTuple):
    readonly : bool = False
    display : Optional[str] = None

MetaFieldOverridesDict = Dict[str, MetaFieldOverrides]

class Field(NamedTuple):
    overrides : MetaFieldOverrides
    name : str
    parse: Callable[[str], Any]
    field_type: Type[Any]

TTuple = TypeVar('TTuple', bound=NamedTuple)

def editable_fields(
    tuple_class : Type[TTuple],
    tuple_fields_overrides : Optional[MetaFieldOverridesDict],
    values : List[Optional[TTuple]],
    default_item: Optional[TTuple],
) -> List[Field]:

    # Rules:
    # values can have None if default_item is not None
    # All values must be an isntance of type_class or None
    for value in values:
        if (value is None and default_item is None) \
            or (value is not None and not isinstance(value, tuple_class)):
            raise ValueError(f"All values in {values} must be an instance of {tuple_class}")

    if default_item is not None and not isinstance(default_item, tuple_class):
        raise ValueError(f"The default item {default_item} must be an instance of {tuple_class}")

    overrides = tuple_fields_overrides or {}

    return [
        Field(
            name = field,
            parse = parser_for(ty),
            overrides = overrides.get(field, MetaFieldOverrides()),
            field_type = ty
        )
        for field, ty in tuple_class.__annotations__.items()
            if can_edit(ty)
    ]

class NamedTupleEditorModel(QAbstractTableModel, Generic[TTuple]):

    def __init__(
        self,
        tuple_class : Type[TTuple],
        tuple_class_overrides : Optional[MetaFieldOverridesDict],
        values : List[Optional[TTuple]],
        panic: Callable[[Exception], None],
        default_item: Optional[TTuple] = None
    ) -> None:
        super().__init__()
        self.__values = values
        self.__fields = editable_fields(
            tuple_class,
            tuple_class_overrides,
            values,
            default_item
        )
        self.__panic = panic
        self.__default_item = default_item

    @property
    def current_values(self):
        return self.__values

    def rowCount(self, parent : Any = None) -> int:
        return len(self.__fields)

    def columnCount(self, parent : Any = None) -> int:
        return len(self.__values)

    def headerData(
            self,
            section : int,
            orientation : Qt.Orientation,
            role: int = Qt.ItemDataRole.DisplayRole
        ):
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return f"{section + 1}"
        elif role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Vertical:
            return self.__get_field_display(section)

        return super().headerData(section, orientation, role)

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        if self.__is_value_checkable(index):
            return super().flags(index) | Qt.ItemFlag.ItemIsUserCheckable

        if self.__is_value_editable(index):
            return super().flags(index) | Qt.ItemFlag.ItemIsEditable

        return super().flags(index)

    def __get_field_display(self, index: int) -> str:
        field = self.__fields[index]

        return field.overrides.display or field.name

    def __get_value(self, index: QModelIndex) -> Optional[Any]:
        field = self.__get_field(index)
        item = self.__get_item(index)

        if item is None:
            return None

        return getattr(item, field)

    def __display_role_data(self, index: QModelIndex):
        return self.__get_value(index)

    def __is_value_editable(self, index: QModelIndex):
        return not self.__get_field_definition(index).overrides.readonly

    def __is_value_checkable(self, index: QModelIndex):
        return issubclass(
            self.__get_field_definition(index).field_type,
            CHECKABLE_TYPES_TUPLE
        )

    def __get_item(self, index: QModelIndex) -> Optional[TTuple]:
        return self.__values[index.column()]

    def __set_item(self, index: QModelIndex, item: Optional[TTuple]):
        self.__values[index.column()] = item

    def __get_field_definition(self, index: QModelIndex) -> Field:
        return self.__fields[index.row()]

    def __get_field(self, index: QModelIndex) -> str:
        return self.__get_field_definition(index).name

    def __get_parser(self, index: QModelIndex) -> Callable[[str], Any]:
        return self.__get_field_definition(index).parse

    def __getitem__(self, ix: int) -> Optional[TTuple]:
        return self.__get_item(self.index(0, ix))
    
    def iterate_values(self) -> Iterable[Optional[TTuple]]:
        return (value for value in self.__values)
    
    def remove(self, *values: TTuple):

        for value in values:
            index = self.__values.index(value)
            if index >= 0:
                self.__values.pop(index)

        self.modelReset.emit()
    
    def append(self, *values: Optional[TTuple]):

        for value in values:
            self.__values.append(value)
        self.modelReset.emit()

    def __setitem__(self, ix: int, value: Optional[TTuple]) -> None:
        self.__set_item(
            self.index(0, ix),
            value
        )

        self.dataChanged.emit(
            self.index(0, ix),
            self.index(self.rowCount() -1, ix)
        )

    def setData(
        self,
        index: QModelIndex,
        value : Any,
        role: int = Qt.ItemDataRole.EditRole
    ) -> bool:

        if role == Qt.ItemDataRole.EditRole:
            return self.__set_edit_data(index, value)
        if role == Qt.ItemDataRole.CheckStateRole:
            return self.__set_checked_sate(index, value)
        return False

    def __set_checked_sate(self, index: QModelIndex, value: Any):
        field = self.__get_field(index)
        item = self.__get_item(index) or self.__default_item

        assert item is not None, "Error! Control flow should not reach this point with item being None"

        new_value = True if value == Qt.CheckState.Checked else False
        self.__set_item(index, item._replace(**{field: new_value}))
        self.dataChanged.emit(index, index)
        return True

    def __set_edit_data(self, index: QModelIndex, value: Any):
        field = self.__get_field(index)
        item = self.__get_item(index) or self.__default_item

        assert item is not None, "Error! Control flow should not reach this point with item being None"
        converter = self.__get_parser(index)

        try:
            self.__set_item(index, item._replace(**{field: converter(value)}))
            self.dataChanged.emit(index, index)
        except Exception as e:
            self.__panic(e)
            return False

        return True

    def write_string(self, data : str):
        for row,line in enumerate(data.splitlines()):
            if row > self.rowCount():
                return
            for col, data in enumerate(line.split('\t')):
                if col >= self.columnCount():
                    break

                self.setData(
                    self.index(row, col),
                    data
                )      

        self.modelReset.emit()

    def __checked_state_data(self, index: QModelIndex) -> Optional[Qt.CheckState]:

        if not self.__is_value_checkable(index):
            return None

        value = self.__get_value(index)

        if value:
            return Qt.CheckState.Checked
        else:
            return Qt.CheckState.Unchecked

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return

        if role == Qt.ItemDataRole.DisplayRole:
            return self.__display_role_data(index)
        elif role == Qt.ItemDataRole.CheckStateRole:
            return self.__checked_state_data(index)

        return None
    
def namedtuple_eidtor(
    table: QTableView,
    *values: Optional[TTuple],
    tuple_type : Optional[Type[TTuple]] = None,
    tuple_field_overrides: Optional[MetaFieldOverridesDict] = None,
    default_item : Optional[TTuple] = None
    ) -> NamedTupleEditorModel[TTuple]:

    def panic(e: Exception):
        show_exception(table, e)

    tuple_type = tuple_type or \
        cast(
            Type[TTuple],
            next(value for value in values if value is not None).__class__
        )
    
    overrides = tuple_field_overrides or getattr(tuple_type, OVERRIDES_FIELD, None)

    editor = NamedTupleEditorModel(
        tuple_type,
        overrides,
        list(values),
        panic,
        default_item
    )
    table.setModel(editor)

    return editor