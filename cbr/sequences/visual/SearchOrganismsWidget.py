from PyQt5.QtCore import QModelIndex, QObject, Qt, pyqtSignal
from PyQt5.QtWidgets import QComboBox, QTableWidgetItem, QWidget
from typing import Any, Dict, List, NamedTuple, Optional, Sequence, Set, Tuple, Union

from ...core.Qt.QtCore import AbstractRecordTableModel
from ...core.Qt.QtWidgets import show_error
from ...core.uRx.core import SubscriptionComposite
from ...core.uRx.dsl import Dsl
from ...core.uRx.qt import QtRx
from ...extra.CbrExtraInteractiveHandler import CbrExtraInteractiveManager, CbrExtraInteractiveHandler 

from ..data import InteractiveInput, InteractiveOutput, SearchArg, SearchArgs, SearchResultRecord
from ..interactive import InteractiveMessage
from .Ui_SearchOrganismsWidget import Ui_SearchOrganismsWidget

TableFields = Dict[str, str]
ExtraFields = Dict[str, TableFields]

class SearchOrganismsState(NamedTuple):
    """Due to the asyncronous nature of the interactive organism search, we need to keep track of the
    current searches that are being processed. This class is used to keep track of the current searches
    and the results that have been received so far.
    """

    current_searches: Set['SearchOrganismsState.Search']
    errors: Optional[List[str]] = None
    results: Optional[List[SearchResultRecord]] = None
    extra_fields: TableFields = {}

    class Search(NamedTuple):
        message_id: int
        search_id: str
        extra_fields: TableFields

        def __hash__(self) -> int:
            return (self.message_id, self.search_id).__hash__()

    @classmethod
    def empty(cls) -> 'SearchOrganismsState':
        return SearchOrganismsState(current_searches=set())
    
    def __on_message(self, message: InteractiveMessage) -> 'SearchOrganismsState':
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
            affected = {
                search
                for search in self.current_searches
                if search.message_id in message.input_uids
            }
            searches = self.current_searches.difference(affected)
            return SearchOrganismsState(
                current_searches=searches,
                errors=[message.error] if len(affected) > 0 else None
            )

        payload = message.payload
        assert payload is not None

        if payload.search_result is not None:
            search_result = payload.search_result
            completed = {
                search
                for search in self.current_searches
                if search.search_id == search_result.search_id
            }
            searches = self.current_searches.difference(completed)

            return SearchOrganismsState(
                current_searches=searches,
                results=payload.search_result.records,
                errors=[error.message for error in search_result.errors],
                extra_fields={
                    k:v
                    for search in completed
                    for k,v in search.extra_fields.items()
                }
            )

        raise ValueError(f"The received message was not understood {message}")

    def is_busy(self) -> bool:
        return len(self.current_searches) > 0

    @classmethod
    def accumulator(
        cls,
        state: 'SearchOrganismsState', 
        value: Union[List['SearchOrganismsState.Search'], InteractiveMessage]
    ) -> 'SearchOrganismsState':
        if isinstance(value, list):
            return SearchOrganismsState(
                current_searches=state.current_searches.union({search for search in value})
            )
        else:
            return state.__on_message(value)

class SaveSearchState(NamedTuple):
    current_saves: Set[int]
    errors: Optional[List[str]] = None

    class Save(NamedTuple):
        message_id: int

    @classmethod
    def empty(cls) -> 'SaveSearchState':
        return SaveSearchState(current_saves=set())

    def is_busy(self) -> bool:
        return len(self.current_saves) > 0

    def __on_message(self, message: InteractiveMessage) -> 'SaveSearchState':

        affected = {
            save_id
            for save_id in self.current_saves
            if save_id in message.input_uids
        }
        errors = None

        if message.error is not None:
            errors = [message.error] if len(affected) > 0 else None
        elif message.payload is not None and message.payload.save_search_result is not None:
            save_search = message.payload.save_search_result
            errors = save_search.errors

        return SaveSearchState(
            current_saves=self.current_saves.difference(affected),
            errors=errors
        )

    @classmethod
    def accumulator(cls, state: 'SaveSearchState', value: Union['SaveSearchState.Save', InteractiveMessage]) -> 'SaveSearchState':
        if isinstance(value, SaveSearchState.Save):
            return SaveSearchState(current_saves=state.current_saves.union({value.message_id}))
        else:
            return state.__on_message(value)


class SearchOrganismsInteractiveHandler(QObject):

    __on_search = pyqtSignal(object)
    __on_save = pyqtSignal(object)

    def __init__(
        self,
        parent: QObject,
        handler: CbrExtraInteractiveHandler[InteractiveOutput, InteractiveInput]
    ):
        super(SearchOrganismsInteractiveHandler, self).__init__(parent)
        self.__current_id = 0
        self.__handler = handler
        self.__search_state = QtRx.observe_signal(self, self.__on_search) \
            .merge_union(handler.observe_message()) \
            .scan(SearchOrganismsState.empty(), SearchOrganismsState.accumulator)

        self.__save_state = QtRx.observe_signal(self, self.__on_save) \
            .merge_union(handler.observe_message()) \
            .scan(SaveSearchState.empty(), SaveSearchState.accumulator)

        self.__subscriptions = SubscriptionComposite([])

    def search_state(self) -> Dsl[SearchOrganismsState]:
        return self.__search_state

    def save_state(self) -> Dsl[SaveSearchState]:
        return self.__save_state

    def busy_state(self) -> Dsl[bool]:
        return self.__search_state.map(lambda s: s.is_busy()).prepend(False) \
            .zip_scan(self.__save_state.map(lambda s: s.is_busy()).prepend(False)) \
            .map(lambda v: v[0] or v[1])

    def search(self, search_args: SearchArgs, extra_fields: ExtraFields) -> None:
        search_task = self.__handler.send_message(InteractiveInput(search=search_args))
        self.__on_search.emit([
            SearchOrganismsState.Search(
                message_id=search_task.message_uid,
                search_id=search.search_id,
                extra_fields=extra_fields.get(search.search_id, {})
            )
            for search in search_args.searches
        ])

    def save(self, records: List[SearchResultRecord]) -> None:
        save_task = self.__handler.send_message(InteractiveInput(save_search=records))
        self.__on_save.emit([
            SaveSearchState.Save(message_id=save_task.message_uid)
        ])

    def next_id(self) -> int:
        self.__current_id += 1
        return self.__current_id

class SearchResultsModelEntry:
    record: SearchResultRecord
    selected: bool
    extra_fields: TableFields

    def __init__(
        self,
        record: SearchResultRecord,
        selected: bool,
        extra_fields: TableFields
    ):
        self.record = record
        self.selected = selected
        self.extra_fields = extra_fields

class SearchResultsModel(AbstractRecordTableModel[str, SearchResultsModelEntry]):
    K_SELECTED = "Selected"
    K_ACCESSION = "Accession"
    K_ORGANISM = "Organism"
    K_TAXID = "TaxId"

    def __init__(self, items: Optional[List[SearchResultsModelEntry]] = None):
        items = items if items is not None else []
        super(SearchResultsModel, self).__init__(items)

        self.__columns = [self.K_SELECTED, self.K_ACCESSION, self.K_ORGANISM, self.K_TAXID]
        self.__extra_fields: List[str] = []

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
        return len(self.__columns) + len(self.__extra_fields)

    def __get_field_for_column(self, column: int) -> Optional[str]:
        if column < len(self.__columns):
            return None
        else:
            ix = column - len(self.__columns)
            return self.__extra_fields[ix] if ix < len(self.__extra_fields) else None

    def headerData(
        self,
        section: int,
        orientation: int,
        role: int = Qt.DisplayRole
    ) -> Any:
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section < len(self.__columns):
                return self.__columns[section]
            else:
                return self.__get_field_for_column(section)

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:

        flags = super().flags(index)
        if self.__columns.index(self.K_SELECTED) == index.column():
            return Qt.ItemIsUserCheckable | flags

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

    def add_results(self, *records: SearchResultRecord, extra_fields: TableFields) -> None:

        new_fields = [field for field in extra_fields if field not in self.__extra_fields]

        if len(new_fields) > 0:
            last = self.columnCount()
            self.beginInsertColumns(QModelIndex(), self.columnCount(), self.columnCount() + len(new_fields) - 1)
            self.__extra_fields += new_fields
            self.endInsertColumns()

        inserts = 0
        for record in records:
            entry = SearchResultsModelEntry(record=record, selected=True, extra_fields=extra_fields)
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

        extra_field = self.__get_field_for_column(column)
        if extra_field is not None:
            return entry.extra_fields.get(extra_field, "")

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

        if column == self.__columns.index(self.K_SELECTED):
            return Qt.Checked if record.selected else Qt.Unchecked
        else:
            return None

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
        manager: CbrExtraInteractiveManager,
        parent: Optional[QWidget] = None
    ):
        super(SearchOrganismsWidget, self).__init__(parent)

        handler: CbrExtraInteractiveHandler[InteractiveOutput, InteractiveInput] = manager.message_handler(
            serializer=InteractiveInput.to_json_dict,
            parser=InteractiveOutput.from_json_dict
        )

        self.setupUi(self)
        self.__handler = SearchOrganismsInteractiveHandler(self, handler)
        self.__subscriptions = SubscriptionComposite([
            QtRx.foreach_signal(self, self.search_button.clicked, self.__on_search_clicked, [bool]),
            QtRx.foreach_signal(self, self.search_table.itemChanged, self.__on_search_table_item_changed, [QTableWidgetItem]),
            self.__handler.search_state().subscribe(QtRx.foreach_ui_observer(self, self.__on_search_state)),
            self.__handler.save_state().subscribe(QtRx.foreach_ui_observer(self, self.__on_save_search_state)),
            self.__handler.busy_state().subscribe(QtRx.foreach_ui_observer(self, self.__set_is_busy))
        ])

    def setupUi(self, SearchOrganismsWidget: QWidget) -> None:
        super().setupUi(SearchOrganismsWidget) #pyright: ignore

        self.busy_progress.setVisible(False)
        self.__results_model = SearchResultsModel()
        self.results_table.setModel(self.__results_model)
        self.busy_progress.setMinimum(0)
        self.busy_progress.setMaximum(0)

    def __set_is_busy(self, is_busy: bool) -> None:
        print("Is busy", is_busy)
        self.busy_progress.setVisible(is_busy)
        self.search_button.setEnabled(not is_busy)
        self.add_to_database_button.setEnabled(not is_busy)

    def __on_save_search_state(self, state: SaveSearchState) -> None:
        errors = state.errors
        if errors is not None and len(errors) > 0:
            raise Exception("\n".join(errors))

    def __on_search_state(self, state: SearchOrganismsState) -> None:

        print("Search state", state)

        results = state.results
        if results is not None:
            self.__results_model.add_results(*results, extra_fields=state.extra_fields)

        errors = state.errors
        if errors is not None and len(errors) > 0:
            raise Exception("\n".join(errors))

    def __on_search_table_item_changed(self, item: Any) -> None:
        if item.row() == 0:
            self.__update_columns_combos()

    def __get_row_fields(self, row: int) -> TableFields:
        search_columns = [
            self.taxid_combo.currentData(),
            self.name_combo.currentData(),
            self.accession_combo.currentData()
        ]

        headers = [
            None if item is None or item.text() == "" else item.text()
            for i in range(self.search_table.columnCount()) 
            for item in [self.search_table.item(0, i)]
        ]

        return {
            header: item.text()
            for i in range(1, self.search_table.columnCount())
                if i not in search_columns
            for header in [headers[i]]
                if header is not None
            for item in [self.search_table.item(row, i)]
                if item is not None
        }

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

    def __get_current_search_args(self) -> List[Tuple[SearchArg, TableFields]]:
        search_terms = [
            (search_arg, self.__get_row_fields(row))
            for row in range(1, self.search_table.rowCount())
            for search_arg in [self.__get_search_arg(row)]
            if search_arg is not None
        ]

        if len(search_terms) == 0:
            raise ValueError("The input does not constitute a valid search. Please ensure that the table contains more than one row and that at least one column is selected.")

        return search_terms

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
        text = item.text() if item is not None else None
        return text if text != "" else None

    def __on_search_clicked(self, _value: Any) -> None:
        args = self.__get_current_search_args()

        search_args = SearchArgs(searches=[arg for arg, _fields in args])
        extra_fields: ExtraFields = {arg.search_id: fields for arg, fields in args}

        self.__handler.search(search_args, extra_fields)

    def set_search_text(self, text: str) -> None:

        rows = text.split("\n")
        self.search_table.setRowCount(len(rows))

        for i, row in enumerate(rows):
            columns = row.split("\t")
            for j, column in enumerate(columns):
                if self.search_table.columnCount() <= j:
                    self.search_table.insertColumn(j)
                self.search_table.setItem(i, j, QTableWidgetItem(column))

