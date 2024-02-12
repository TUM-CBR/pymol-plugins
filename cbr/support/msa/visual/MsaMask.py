from Bio.Align import MultipleSeqAlignment
from enum import Enum
from typing import NamedTuple, List, Set, Tuple

from ...display.sequence import RESIDUE_COLORS

class MaskPositionMode(Enum):
    HIGHLIGHT = 0
    HIDE = 1

class MsaMask(NamedTuple):
    masked_rows : Set[int]
    masked_columns : Set[int]
    row_mappings : List[int]
    col_mappings : List[int]
    mask_mode : MaskPositionMode

    @property
    def row_count(self):
        return len(self.row_mappings)

    @property
    def col_count(self):
        return len(self.col_mappings)

    def get_position(self, row : int, col : int) -> Tuple[int, int]:
        return (
            self.row_mappings[row],
            self.col_mappings[col]
        )

    def is_masked(self, row : int, col : int) -> bool:
        return row in self.masked_rows or col in self.masked_columns

    @staticmethod
    def default_mask(alignment: MultipleSeqAlignment) -> 'MsaMask':

        return MsaMask(
            masked_columns=set([]),
            masked_rows=set([]),
            row_mappings=[i for i in range(0, len(alignment))],
            col_mappings=[i for i in range(0, alignment.get_alignment_length())],
            mask_mode=MaskPositionMode.HIDE
        )