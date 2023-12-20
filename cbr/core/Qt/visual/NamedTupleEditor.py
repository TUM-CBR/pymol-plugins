from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from PyQt5.QtWidgets import QTableView
from typing import Any, Callable, cast, List, NamedTuple, Optional, Type

from ..QtWidgets import show_exception

ParseValue = Callable[[str], Any]

EDITABLE_TYPES = {
    str: lambda s: s,
    float: lambda s: float(s),
    int: lambda s: int(s)
}

EDITABLE_TYPES_TUPLE = tuple(EDITABLE_TYPES.keys())

def can_edit(v: Any) -> bool:

    if not isinstance(v, Type):
        v = v.__class__

    return issubclass(v, EDITABLE_TYPES_TUPLE)

def parser_for(value: Type) -> ParseValue:

    for ty, parser in EDITABLE_TYPES.items():
        if issubclass(value, ty):
            return parser

    raise Exception(f"There is no editor for value {value}")

OVERRIDES_FIELD = '__editor_field_overrides__'

class MetaFieldOverrides(NamedTuple):
    readonly : bool = False
    display : Optional[str] = None

class Field(NamedTuple):
    overrides : MetaFieldOverrides
    name : str
    parse: Callable[[str], Any]

def editable_fields(
    tuple_class : Type[NamedTuple],
    values : List[Optional[NamedTuple]],
    default_item: Optional[NamedTuple]
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

    overrides = getattr(tuple_class, OVERRIDES_FIELD, {})

    return [
        Field(
            name = field,
            parse=parser_for(ty),
            overrides=overrides.get(field, MetaFieldOverrides())
        )
        for field, ty in tuple_class.__annotations__.items()
            if can_edit(ty)
    ]

class NamedTupleEditorModel(QAbstractTableModel):

    def __init__(
        self,
        tuple_class : Type[NamedTuple],
        values : List[Optional[NamedTuple]],
        panic: Callable[[Exception], None],
        default_item: Optional[NamedTuple] = None
    ) -> None:
        super().__init__()
        self.__values = values
        self.__fields = editable_fields(tuple_class, values, default_item)
        self.__panic = panic
        self.__default_item = default_item

    @property
    def current_values(self):
        return self.__values

    def rowCount(self, parent = None) -> int:
        return len(self.__fields)

    def columnCount(self, parent = None) -> int:
        return len(self.__values)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return f"{section + 1}"
        elif role == Qt.DisplayRole and orientation == Qt.Vertical:
            return self.__get_field_display(section)

        return super().headerData(section, orientation, role)

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        if self.__is_value_editable(index):
            return super().flags(index) | Qt.ItemIsEditable

        return super().flags(index)

    def __get_field_display(self, index: int) -> str:
        field = self.__fields[index]

        return field.overrides.display or field.name

    def __display_role_data(self, index: QModelIndex):
        field = self.__get_field(index)
        item = self.__get_item(index)

        return item and getattr(item, field)

    def __is_value_editable(self, index: QModelIndex):
        return not self.__get_field_definition(index).overrides.readonly

    def __get_item(self, index: QModelIndex) -> Optional[NamedTuple]:
        return self.__values[index.column()]

    def __set_item(self, index: QModelIndex, item: Optional[NamedTuple]):
        self.__values[index.column()] = item

    def __get_field_definition(self, index: QModelIndex) -> Field:
        return self.__fields[index.row()]

    def __get_field(self, index: QModelIndex) -> str:
        return self.__get_field_definition(index).name

    def __get_parser(self, index: QModelIndex) -> Callable[[str], Any]:
        return self.__get_field_definition(index).parse

    def setData(self, index: QModelIndex, value, role = Qt.EditRole) -> bool:

        if role != Qt.EditRole:
            return False

        field = self.__get_field(index)
        item = self.__get_item(index) or self.__default_item

        assert item is not None, "Error! Control flow should not reach this point with item being None"
        current = getattr(item, field)
        converter = self.__get_parser(index)

        try:
            self.__set_item(index, item._replace(**{field: converter(value)}))
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

    def data(self, index: QModelIndex, role = Qt.DisplayRole):

        if not index.isValid():
            return

        if role == Qt.DisplayRole:
            return self.__display_role_data(index)

        return None
    
def namedtuple_eidtor(
    table: QTableView,
    *values: Optional[NamedTuple],
    tuple_type : Optional[Type[NamedTuple]] = None,
    default_item : Optional[NamedTuple] = None
    ) -> NamedTupleEditorModel:

    def panic(e: Exception):
        show_exception(table, e)

    tuple_type = tuple_type or \
        cast(
            Type[NamedTuple],
            next(value for value in values if value is not None).__class__
        )

    editor = NamedTupleEditorModel(
        tuple_type,
        list(values),
        panic,
        default_item
    )
    table.setModel(editor)

    return editor