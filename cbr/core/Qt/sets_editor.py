from PyQt5.QtCore import QAbstractTableModel, QModelIndex, Qt
from typing import Any, Callable, Generic, List, Optional, Set, TypeVar

TSetItem = TypeVar('TSetItem')

class SetsEditorModel(QAbstractTableModel, Generic[TSetItem]):

    def __init__(
        self,
        universe: Set[TSetItem],
        sets: List[Set[TSetItem]],
        exclusive: bool = False,
        render: Optional[Callable[[TSetItem], str]] = None
    ):
        super().__init__()
        self.__sets = sets
        self.__exclusive = exclusive
        self.__render_fn = render if render is not None else str
        self.__set_universe(universe)

    def __set_universe(self, universe: Set[TSetItem]):
        self.__universe = list(universe)
        self.__sets = [
            set(s for s in selected if s in universe)
            for selected in self.__sets
        ]
        self.__exclusive_assignments: List[Optional[int]] = [
            next(
                (i for i,set in enumerate(self.__sets) if item in set),
                None
            )
            for item in self.__universe
        ]

        self.modelReset.emit()

    def __render_text(self, item: TSetItem) -> str:
        return self.__render_fn(item)

    def rowCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__universe)
    
    def columnCount(self, parent: Optional[QModelIndex] = None) -> int:
        return len(self.__sets)
    
    def headerData(
        self,
        section: int,
        orientation: Qt.Orientation,
        role: int = Qt.ItemDataRole.DisplayRole
    ) -> Any:
        
        if role != Qt.ItemDataRole.DisplayRole \
            or orientation == Qt.Orientation.Horizontal:
            return super().headerData(section, orientation, role)
        
        return self.__render_text(self.__universe[section])
    
    def __get_universe_item_index(self, index: QModelIndex) -> int:
        return index.row()
    
    def __get_set_index(self, index: QModelIndex) -> int:
        return index.column()
    
    def __get_set_for_index(self, index: QModelIndex) -> Set[TSetItem]:
        return self.__sets[self.__get_set_index(index)]
    
    def __get_universe_item(self, index: QModelIndex) -> TSetItem:
        return self.__universe[self.__get_universe_item_index(index)]
    
    def __set_exclusive(self, index: QModelIndex) -> None:
        set_ix = self.__get_set_index(index)
        item_ix = self.__get_universe_item_index(index)
        self.__exclusive_assignments[item_ix] = set_ix

    def __can_be_assigned(self, index: QModelIndex) -> bool:
        """ Determine if the index is a valid assignment. When
        the "exclusive" option is enabled, every value in the
        universe can belong at most one to a set.
        """

        if not self.__exclusive:
            return True

        member_ix = self.__get_universe_item_index(index)
        set_ix = self.__get_set_index(index)
        return self.__exclusive_assignments[member_ix] == set_ix
    
    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index) | Qt.ItemFlag.ItemIsUserCheckable

        if not self.__can_be_assigned(index):
            flags &= ~Qt.ItemFlag.ItemIsEnabled

        return  flags

    def setData(self, index: QModelIndex, value: Any, role: int = Qt.ItemDataRole.EditRole) -> bool:

        if not index.isValid() or role != Qt.ItemDataRole.CheckStateRole:
            return False
        
        set_to_update = self.__get_set_for_index(index)
        item_to_add = self.__get_universe_item(index)
        set_to_update.add(item_to_add)

        if self.__exclusive:
            # Essentially, all columns but one will need an update
            # so we may as well just reset the whole thing
            self.__set_exclusive(index)
            self.modelReset.emit()
        return super().setData(index, value, role)
    
    def __checkd_state_data(self, index: QModelIndex) -> Qt.CheckState:
        index_set = self.__get_set_for_index(index)
        item = self.__get_universe_item(index)

        if item in index_set:
            return Qt.CheckState.Checked
        else:
            return Qt.CheckState.Unchecked
    
    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole) -> Any:

        if not index.isValid():
            return None
        
        if role == Qt.ItemDataRole.CheckStateRole:
            return self.__checkd_state_data(index)
        return super().data(index, role)