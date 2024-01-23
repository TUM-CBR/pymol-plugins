from abc import abstractmethod
from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from typing import Any, Callable, Dict, Optional, Type

ParseValue = Callable[[str], Any]

EDITABLE_TYPES : Dict[Type[Any], Callable[[str], Any]] = {
    str: lambda s: s,
    float: lambda s: float(s),
    int: lambda s: int(s),
    bool: lambda s: bool(s)
}

CHECKABLE_TYPES_TUPLE = (bool,)

EDITABLE_TYPES_TUPLE = tuple(EDITABLE_TYPES.keys())

def can_edit(v: Any) -> bool:

    if not isinstance(v, Type):
        v = v.__class__

    return issubclass(v, EDITABLE_TYPES_TUPLE)

def parser_for(value: Type[Any]) -> ParseValue:

    for ty, parser in EDITABLE_TYPES.items():
        if issubclass(value, ty):
            return parser

    raise Exception(f"There is no editor for value {value}")

class EditRecordModelBase(QAbstractTableModel):

    @abstractmethod
    def __get_value__(self, index: QModelIndex) -> Optional[Any]:
        raise NotImplementedError("__get_value__ needs to be overriden")
    
    @abstractmethod
    def __set_value__(self, index: QModelIndex, value: Any) -> bool:
        raise NotImplementedError("__set_value__ needs to be overriden")
    
    def __get_type__(self, index: QModelIndex) -> Optional[Type[Any]]:
        value = self.__get_value__(index)

        if value is None:
            return None
        
        return type(value)

    def __set_value(self, index: QModelIndex, value: Any) -> bool:

        value_type = self.__get_type__(index)

        if value_type is None:
            return False
        
        parser = parser_for(value_type)
        return self.__set_value__(index, parser(value))
    
    def __can_edit__(self, index: QModelIndex):
        return can_edit(self.__get_type__(index))
    
    def setData(
        self,
        index : QModelIndex,
        value : Any,
        role : int = Qt.ItemDataRole.EditRole
    ) -> bool:
        
        if role == Qt.ItemDataRole.EditRole:
            result = self.__set_value(index, value)
            self.dataChanged.emit(index, index)
            return result
        
        return False
