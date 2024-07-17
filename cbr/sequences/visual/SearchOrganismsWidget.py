from PyQt5.QtCore import QModelIndex, QObject, Qt, pyqtSignal
from PyQt5.QtWidgets import QComboBox, QWidget
from typing import Any, List, NamedTuple, Optional, Set, Union

from ...core.Qt.QtCore import AbstractRecordTableModel
from ...core.Qt.QtWidgets import show_error
from ...core.uRx.core import SubscriptionComposite
from ...core.uRx.dsl import Dsl
from ...core.uRx.qt import QtRx
from ...extra.CbrExtraInteractiveHandler import CbrExtraInteractiveHandler, CbrMessage 

from ..data import InteractiveInput, InteractiveOutput, SearchArg, SearchArgs, SearchResultRecord
from .Ui_SearchOrganismsWidget import Ui_SearchOrganismsWidget

class OrganismSearchState(NamedTuple):
    message_id: int
    search_id: str

SearchResultMessage = CbrMessage[InteractiveOutput, InteractiveInput]

class SearchOrganismsState(NamedTuple):
    """Due to the asyncronous nature of the interactive organism search, we need to keep track of the
    current searches that are being processed. This class is used to keep track of the current searches
    and the results that have been received so far.
    """

    current_searches: Set[OrganismSearchState]
    errors: Optional[List[str]] = None
    results: Optional[List[SearchResultRecord]] = None 

    @classmethod
    def empty(cls) -> 'SearchOrganismsState':
        return SearchOrganismsState(current_searches=set())
    
    def on_message(self, message: SearchResultMessage) -> 'SearchOrganismsState':
        """This method is called every time a search result message is received. It
        is responsible for updating the state of the search based on the message that
        was received.

        The possible outcomes are:
        - The message contains an error, in which case the error is added to the state
            and we also remove the queued searches matching the message uid of the error message.
        - The message contains search results. In that case, we remove the search from the queue
            and set the results property of the state to contain the latest results. We also add
            any errors that occurred during the search.
        """

        if message.error is not None:
            searches = self.current_searches.difference({
                search
                for search in self.current_searches
                if search.message_id in message.input_uids
            })
            return SearchOrganismsState(
                current_searches=searches,
                errors=[message.error]
            )

        payload = message.payload
        assert payload is not None

        if payload.search_result is not None:
            search_result = payload.search_result
            searches = self.current_searches.difference({
                search
                for search in self.current_searches
                if search.search_id == search_result.search_id
            })

            return SearchOrganismsState(
                current_searches=searches,
                results=payload.search_result.records,
                errors=[error.message for error in search_result.errors]
            )

        raise ValueError(f"The received message was not understood {message}")

    def is_done(self) -> bool:
        return len(self.current_searches) == 0

    @classmethod
    def accumulator(
        cls,
        state: 'SearchOrganismsState', 
        value: Union[OrganismSearchState, SearchResultMessage]
    ) -> 'SearchOrganismsState':
        if isinstance(value, OrganismSearchState):
            return SearchOrganismsState(current_searches=state.current_searches.union({value}))
        else:
            return state

class SearchOrganismsInteractiveHandler(QObject):

    __on_search = pyqtSignal(List[OrganismSearchState])

    def __init__(
        self,
        parent: QObject,
        handler: CbrExtraInteractiveHandler[InteractiveOutput, InteractiveInput],
    ):
        super(SearchOrganismsInteractiveHandler, self).__init__(parent)
        self.__current_id = 0
        self.__handler = handler

        self.__search_state = QtRx.observe_signal(self, self.__on_search) \
            .merge_union(handler.observe_message()) \
            .scan(SearchOrganismsState.empty(), SearchOrganismsState.accumulator)

        self.__subscriptions = SubscriptionComposite([])

    def search_state(self) -> Dsl[SearchOrganismsState]:
        return self.__search_state

    def search(self, search_args: SearchArgs) -> None:
        search_task = self.__handler.send_message(InteractiveInput(search=search_args))
        self.__on_search.emit([
            OrganismSearchState(message_id=search_task.message_uid, search_id=search.search_id)
            for search in search_args.searches
        ])

    def next_id(self) -> int:
        self.__current_id += 1
        return self.__current_id

class SearchResultsModelEntry:
    record: SearchResultRecord
    selected: bool

    def __init__(
        self,
        record: SearchResultRecord,
        selected: bool
    ):
        self.record = record
        self.selected = selected

class SearchResultsModel(AbstractRecordTableModel[str, SearchResultsModelEntry]):
    K_SELECTED = "Selected"
    K_ACCESSION = "Accession"
    K_ORGANISM = "Organism"
    K_TAXID = "TaxId"

    def __init__(self, items: Optional[List[SearchResultsModelEntry]] = None):
        items = items if items is not None else []
        super(SearchResultsModel, self).__init__(items)

        self.__columns = [self.K_SELECTED, self.K_ACCESSION, self.K_ORGANISM, self.K_TAXID]

    def get_uid(self, value: SearchResultsModelEntry) -> str:
        return value.record.accession

    def rowCount(
        self,
        parent: Optional[QModelIndex] = None #pyright: ignore
    ) -> int:
        return self.get_record_count()

    def columnCount(
        self,
        parent: Optional[QModelIndex] = None #pyright: ignore
    ) -> int:
        return len(self.__columns)

    def headerData(
        self,
        section: int,
        orientation: int,
        role: int = Qt.DisplayRole
    ) -> Any:
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            return self.__columns[section]

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)
        if self.__columns.index(self.K_SELECTED) == index.column():
            return Qt.ItemIsUserCheckable | Qt.ItemIsEnabled | flags

        return flags

    def setData(
        self,
        index: QModelIndex,
        value: Any,
        role: int = Qt.EditRole
    ) -> bool:

        if role == Qt.CheckStateRole:
            self.get_record(index.row()).selected = value == Qt.Checked
            self.dataChanged.emit(index, index)
            return True

        return False

    def add_results(self, *records: SearchResultRecord) -> None:

        inserts = 0
        for record in records:
            entry = SearchResultsModelEntry(record=record, selected=True)
            if not self.has_record(entry):
                self.add_records(entry)
                inserts += 1

        if inserts > 0:
            self.beginInsertRows(QModelIndex(), self.get_record_count(), self.get_record_count() + inserts - 1)
            self.endInsertRows()

    def __data_display_role(
        self,
        entry: SearchResultsModelEntry,
        column: int
    ) -> Any:
        name = self.__columns[column]
        record = entry.record

        if name == self.K_ACCESSION:
            return record.accession
        elif name == self.K_ORGANISM:
            return record.organism.name
        elif name == self.K_TAXID:
            return record.organism.taxid
        else:
            return None

    def __data_check_state_role(
        self,
        record: SearchResultsModelEntry,
        column: int
    ) -> Any:
        return Qt.Checked if record.selected else Qt.Unchecked

    def data(
        self,
        index: QModelIndex,
        role: int = Qt.DisplayRole
    ) -> Any:
    
        if not index.isValid():
            return None
    
        row = index.row()
        column = index.column()
        record = self.get_record(row)

        if role == Qt.DisplayRole:
            return self.__data_display_role(record, column)
        elif role == Qt.CheckStateRole:
            return self.__data_check_state_role(record, column)
        else:
            return None 

class SearchOrganismsWidget(QWidget, Ui_SearchOrganismsWidget):
    def __init__(
        self,
        handler: CbrExtraInteractiveHandler[InteractiveOutput, InteractiveInput],
        parent: Optional[QWidget] = None
    ):
        super(SearchOrganismsWidget, self).__init__(parent)
        self.setupUi(self)
        self.__handler = SearchOrganismsInteractiveHandler(self, handler)
        self.__subscriptions = SubscriptionComposite([
            QtRx.foreach_signal(self, self.search_button.clicked, self.__on_search_clicked),
            QtRx.foreach_signal(self, self.search_table.itemChanged, self.__on_search_table_item_changed),
            self.__handler.search_state().for_each(self.__on_search_state)
        ])

    def setupUi(self, SearchOrganismsWidget: QWidget) -> None:
        super().setupUi(SearchOrganismsWidget) #pyright: ignore

        self.busy_progress.setVisible(False)
        self.__results_model = SearchResultsModel()
        self.results_table.setModel(self.__results_model)

    def __set_is_searching(self, is_searching: bool) -> None:
        self.busy_progress.setVisible(is_searching)

    def __on_search_state(self, state: SearchOrganismsState) -> None:
        self.__set_is_searching(not state.is_done())

        results = state.results
        if results is not None:
            self.__results_model.add_results(*results)

        errors = state.errors
        if errors is not None:
            show_error(
                self,
                "Search Errors",
                "\n".join(errors)
            )

    def __on_search_table_item_changed(self, item: Any) -> None:
        if item.row() == 0:
            self.__update_columns_combos()

    def __get_search_arg(self, row: int) -> Optional[SearchArg]:
        taxid = self.__get_combo_item(self.taxid_combo, row)
        name = self.__get_combo_item(self.name_combo, row)
        accession = self.__get_combo_item(self.accession_combo, row)

        if taxid is None and name is None and accession is None:
            return None

        search_id = self.__handler.next_id()
        return SearchArg(
            search_id = f"search_{search_id}",
            tax_ids = [int(taxid)] if taxid is not None else None,
            names = [name] if name is not None else None,
            accession = accession
        )

    def __get_current_search_args(self) -> SearchArgs:
        return SearchArgs([
            search_arg
            for row in range(self.search_table.rowCount())
            for search_arg in [self.__get_search_arg(row)]
            if search_arg is not None
        ])

    def __update_columns_combos(self) -> None:

        columns = [
            (i, item.text())
            for i in range(self.search_table.columnCount())
            for item in [self.search_table.item(0, i)]
                if item is not None and item.text() != ""
        ]

        combos = [self.name_combo, self.accession_combo, self.taxid_combo]
        for combo in combos:
            combo.clear()
            combo.addItem("<none>", None)

            for i, name in columns:
                combo.addItem(name, i)

    def __get_combo_item(self, combo: QComboBox, row: int) -> Optional[str]:
        selection = combo.currentData()

        if selection is None:
            return None

        item = self.search_table.item(row, selection)
        return item.text() if item is not None else None

    def __on_search_clicked(self, _value: Any) -> None:
        args = self.__get_current_search_args()

        self.__handler.search(args)
