from typing import Callable, Optional
import pymol
from PyQt5.QtCore import QAbstractItemModel, Qt
from PyQt5.QtCore import QItemSelection, QItemSelectionModel, pyqtSlot
from PyQt5.QtWidgets import QTableView

from ....control import viter

class PymolResidueSelectionModel(QItemSelectionModel):

    def __init__(
        self,
        view : QTableView,
        get_selection_scope: Callable[[], Optional[str]],
        resi_role : Qt.ItemDataRole,
        selection_name : str = "sele"
    ):
        super().__init__()

        self.__selection_name = selection_name
        self.__role = resi_role
        self.__get_selection_scope = get_selection_scope
        self.__view = view
        self.selectionChanged.connect(self.__on_selection_changed)

    def setModel(self, model: QAbstractItemModel):
        super().setModel(model)
        self.__view.setModel(model)
        self.__view.setSelectionModel(self)

    def __clear_selection(self):
        selections = pymol.cmd.get_names("selections")
        if self.__selection_name in selections:
            pymol.cmd.select(self.__selection_name)

    @pyqtSlot(QItemSelection,QItemSelection)
    def __on_selection_changed(self, new: QItemSelection, old: QItemSelection):

        if self.__get_selection_scope() is None:
            return

        self.__clear_selection()

        resi = [
            f"resi {position}"
            for item in self.selectedIndexes()
            for position in viter(item.data(self.__role))
        ]

        if len(resi) == 0:
            return

        selection = " or ".join(resi)
        pymol.cmd.select(
            self.__selection_name,
            f"{self.__get_selection_scope()} and ({selection})"
        )

pymol_residue_selection_model = PymolResidueSelectionModel