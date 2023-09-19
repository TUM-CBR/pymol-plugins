import itertools
from pymol import cmd
from PyQt5.QtCore import QModelIndex, QObject, pyqtSlot
from PyQt5.QtWidgets import QTableWidget
from typing import Callable, cast, Dict, Iterable, List, NamedTuple, Optional, Set

import pymol

from ... import color

class CellSelectionContext(NamedTuple):
    residues : Set[int]
    context : object

class PymolResidueSelectorBase(QObject):
    """ Helper QObject that can be used in a QTableWidget. It will perform a pymol selection
    whenever a cell in the table is selected. A mapping from table row to residue number is
    provided and that mapping will be used to make the selection.
    """

    def __init__(
        self,
        qtable : QTableWidget,
        scope_selection : Optional[str] = None,
    ):
        super(PymolResidueSelectorBase, self).__init__()
        self.__selection : Optional[str] = scope_selection
        self.__qtable = qtable
        qtable.itemSelectionChanged.connect(self.__on_item_selection_changed)

    def __cell_to_selection_mapping__(
        self,
        table: QTableWidget,
        index: QModelIndex
    ) -> Optional[CellSelectionContext]:
        return None

    def reset(self):
        self.__selection = None
        self.__on_reset__()

    def __on_reset__(self):
        pass

    def set_selection(self, selection : Optional[str]=None):
        self.__selection = selection

    def __on_selection__(self, selection_scope : str, cell_selection : Iterable[CellSelectionContext]):
        pass

    @pyqtSlot()
    def __on_item_selection_changed(self):
        if self.__selection is None:
            return

        selection = self.__qtable.selectedIndexes()
        residues = [
            context
            for index in selection
            for context in [self.__cell_to_selection_mapping__(self.__qtable, index)] if context is not None
        ]

        self.__on_selection__(self.__selection, residues)


class PymolResidueSelector(PymolResidueSelectorBase):
    __counter = itertools.count()
    
    def __init__(
        self,
        qtable : QTableWidget,
        scope_selection : Optional[str] = None,
        selection_name : Optional[str] = None
    ):
        super(PymolResidueSelector, self).__init__(qtable, scope_selection=scope_selection)
        self.__default_selection_name = "TableSelection_%i" % next(PymolResidueSelector.__counter)
        self.__mappings : Dict[int, List[int]] = {}
        self.__selection_name = selection_name
    
    def __cell_to_selection_mapping__(self, table: QTableWidget, index : QModelIndex) -> Optional[CellSelectionContext]:
        mappings = self.__mappings.get(index.row())

        if mappings is not None:
            return CellSelectionContext(residues=set(mappings), context=None)

    def __on_reset__(self):
        super().__on_reset__()
        self.__mappings = {}

    def set_item_riesidues(self, item : int, residue : List[int]):
        self.__mappings[item] = residue

    def __on_selection__(self, selection_scope : str, cell_selection : Iterable[CellSelectionContext]):
        residues_selection = [
            "resi %i" % residue
            for item in cell_selection
            for residue in item.residues
        ]

        selection = self.__selection_name or self.__default_selection_name

        if(len(residues_selection) == 0):
            cmd.select(
                selection
            )
        else:
            selection_str = " & ".join([
                selection_scope,
                " | ".join(residues_selection)
            ])
            cmd.select(
                selection,
                selection=selection_str
            )

class ResidueWithIntensity(NamedTuple):
    residues : Set[int]
    intensity : float

SelectionToResideAndIntensity = Callable[[QTableWidget, QModelIndex], Optional[ResidueWithIntensity]]

class PymolResidueColorSelector(PymolResidueSelectorBase):

    def __init__(
        self,
        qtable : QTableWidget,
        residue_mapping : SelectionToResideAndIntensity,
        scope_selection : Optional[str] = None
    ):

        super(PymolResidueColorSelector, self).__init__(qtable, scope_selection)
        self.__residue_mapping = residue_mapping

    def __cell_to_selection_mapping__(self, table: QTableWidget, index: QModelIndex) -> Optional[CellSelectionContext]:
        selection = self.__residue_mapping(table, index)

        if selection is None:
            return None

        return CellSelectionContext(residues=selection.residues, context=selection.intensity)

    def __on_selection__(self, selection_scope : str, cell_selection : Iterable[CellSelectionContext]):

        structure_colors : Dict[int, int] = {}
        def __count_color__(c: int):
            if c in structure_colors:
                structure_colors[c] = structure_colors[c] + 1
            else:
                structure_colors[c] = 1
        pymol.cmd.iterate(selection_scope, 'count_color(color)', space={'count_color': __count_color__})
        consensus_color = \
            color.pymol_color_tuple_to_pymol_color(
                pymol.cmd.get_color_tuple(
                    max(structure_colors, key=lambda k:structure_colors[k]
                )
            )
        )
        pymol.cmd.color(consensus_color, selection_scope)
        pymol.cmd.label(f"{selection_scope} and name CA", "")

        for i,selection in enumerate(cell_selection):

            resi_color = color.to_pymol_color(
                color.change_intensity(
                    color.get_color_rgb(i),
                    cast(float, selection.context)
                )
            )

            resi_color = color.to_pymol_color(
                    color.get_color_rgb(i)
            )

            for residue in selection.residues:
                region_selection = f"{selection_scope} and (resi {residue})"
                pymol.cmd.color(
                    resi_color,
                    region_selection
                )
                pymol.cmd.label(f"{region_selection} and name CA", str(selection.context))

def pymol_residues_color_selector(
    qtable : QTableWidget,
    residue_mapping : SelectionToResideAndIntensity,
    scope_selection : Optional[str] = None
):
    return PymolResidueColorSelector(qtable, residue_mapping, scope_selection)