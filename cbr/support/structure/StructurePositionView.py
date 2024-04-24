from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import QObject
from PyQt5.QtGui import QColor
from typing import Callable, Dict, Iterable, List, NamedTuple, Optional, Sequence

from ...clustal.Clustal import Clustal
from ...core.pymol.structure import StructureSelection

from .AbstractCompositeTableModel import *
from .StructuresAlignmentMapper import StructureAlignmentEntry, StructuresAlignmentMapper

class SetStructureArg(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str

class PositionsAndColor(NamedTuple):
    positions: Sequence[int]
    structure_color: Optional[QColor] = None

    def map_to_structure(self, msa_mapping: Dict[int, int]) -> Optional['PositionsAndColor']:

        mapped_positions = [msa_mapping[p] for p in self.positions if p in msa_mapping]

        if len(mapped_positions) == 0:
            return None
        else:
            return self._replace(positions=mapped_positions)

class PositionEntry(NamedTuple):
    positions_and_colors: Sequence[PositionsAndColor]

    @property
    def positions(self) -> Sequence[int]:
        return [
            position
            for entry in self.positions_and_colors
            for position in entry.positions
        ]
    
    def map_positions_and_colors(self, msa_mapping: Dict[int, int]) -> Iterable[PositionsAndColor]:

        return (
            mapped
            for p in self.positions_and_colors
            for mapped in [p.map_to_structure(msa_mapping)]
            if mapped is not None
        )
    
    def as_view_attribute(self, msa_mapping: Optional[Dict[int, int]] = None) -> ViewRecordAttributes:

        if msa_mapping is None:
            records = self.positions_and_colors
        else:
            records = self.map_positions_and_colors(msa_mapping)

        display = [
            str(p)
            for mapped in records
            for p in mapped.positions
        ]

        if len(display) == 0:
            return ViewRecordAttributes.default()

        return  ViewRecordAttributes(
            display=",".join(display)
        )
    
    def as_pymol_attribute(
        self,
        structure: StructureSelection,
        msa_mapping: Dict[int, int]
    ) -> PymolRecordAttributes:

        selections: List[StructureSelection] = []
        colors: List[Optional[QColor]] = []

        for mapped in self.map_positions_and_colors(msa_mapping):
            selections.append(structure.scoped(mapped.positions))
            colors.append(mapped.structure_color) # type: ignore

        assert len(selections) == len(colors), "Algorithm is inconsistent"

        if len(selections) == 0:
            return PymolRecordAttributes.default()
        else:
            return PymolRecordAttributes(
                selections=selections,
                colors=colors
            )

class AbstractStructurePositionView(AbstractRecordView[TModelRecord]):
    """
    Table view which provides columns that can map positions in a multiple sequence
    alignment to positions in pymol structures. This view also allows defining a coloring
    scheme that will be applied to the positions in the pymol structure.

    This view will add the following columns:
        "Position (msa)": The position in the multiple squence alignment
        "Position (strucute)": One column per structure with the position of the structure
        matching the position of the alignment.
    """

    __mapper: Optional[StructuresAlignmentMapper]
    __structures: Sequence[StructureAlignmentEntry]
    __msa_to_resv: Sequence[Dict[int, int]]

    K_POSITION = "Position (msa)"

    def __init__(
        self,
        clustal: Clustal,
        parent: Optional[QObject] = None,
        msa: Optional[MultipleSeqAlignment] = None
    ) -> None:
        super().__init__(parent)

        self.__clustal = clustal
        self.__mapper = None
        self.__structures = []
        self.__msa_to_resv = []

        if msa is not None:
            self.set_alignment(msa)

    def __headers(self) -> Sequence[ViewHeaderSpec]:
        structure_headers = [
            f"Position ({structure.structure.show()})"
            for structure in self.__structures
        ]

        return [
            ViewHeaderSpec(name)
            for name in [self.K_POSITION] + structure_headers
        ]

    def set_structures(self, structures: Sequence[SetStructureArg]):

        mapper = self.__mapper

        if mapper is None:
            raise ValueError("Alignemnt must be set before setting the structures")

        self.__structures = [
            mapper.map_structure(
                structure.structure,
                structure.reference_sequence_id
            )
            for structure in structures
        ]

        self.__msa_to_resv = [
            mapping.msa_to_resv()
            for mapping in self.__structures
        ]

        self.content_changed.emit(self)

    def set_alignment(self, msa: MultipleSeqAlignment):

        self.__mapper = StructuresAlignmentMapper(
            msa,
            self.__clustal
        )

        self.set_structures([
            SetStructureArg(
                structure=entry.structure,
                reference_sequence_id=entry.reference_sequence_id
            )
            for entry in self.__structures
        ])

    def __record_to_position__(self, msa: MultipleSeqAlignment, record: TModelRecord) -> Optional[PositionEntry]:
        raise NotImplementedError("The function __record_to_position__ must be implemented.")
    
    def __get_record_attributes(self, mapper: StructuresAlignmentMapper, record: TModelRecord) -> Sequence[ViewRecordAttributes]:

        position = self.__record_to_position__(mapper.msa, record)

        if position is None:
            position_view = self.DEFAULT_QT_ATTRIBUTES
            structures_view = [self.DEFAULT_QT_ATTRIBUTES for _ in self.__msa_to_resv]
            return [position_view] + structures_view

        position_view = position.as_view_attribute()

        structures_view = [
            position.as_view_attribute(mapping)
            for mapping in self.__msa_to_resv
        ]

        return [position_view] + structures_view
    
    def __get_pymol_attributes(self, mapper: StructuresAlignmentMapper, record: TModelRecord) -> Sequence[PymolRecordAttributes]:

        position = self.__record_to_position__(mapper.msa, record)

        if position is None:
            return [self.DEFAULT_PYMOL_ATTRIBUTES] \
                + [self.DEFAULT_PYMOL_ATTRIBUTES for _ in self.__msa_to_resv]

        return [self.DEFAULT_PYMOL_ATTRIBUTES] + \
            [
                position.as_pymol_attribute(structure.structure, mapping)
                for structure, mapping in zip(self.__structures, self.__msa_to_resv)
            ]

    def attributes(self, records: Sequence[TModelRecord]) -> ViewRecords:
        
        mapper = self.__mapper

        if mapper is None:
            return super().attributes(records)
        
        qt_attributes: Sequence[Sequence[ViewRecordAttributes]] = []
        pymol_attributes: Sequence[Sequence[PymolRecordAttributes]] = []

        for record in records:
            qt_attributes.append(self.__get_record_attributes(mapper, record))
            pymol_attributes.append(self.__get_pymol_attributes(mapper, record))

        return ViewRecords(
            headers=self.__headers(),
            qt_attributes=qt_attributes,
            pymol_attributes=pymol_attributes
        )

GetMsaPosition = Callable[[MultipleSeqAlignment, TModelRecord], Optional[PositionEntry]]

class StructurePositionView(AbstractStructurePositionView[TModelRecord]):

    def __init__(
        self,
        clustal: Clustal,
        get_msa_position: GetMsaPosition[TModelRecord],
        parent: Optional[QObject] = None,
        msa: Optional[MultipleSeqAlignment] = None
    ) -> None:
        super().__init__(clustal, parent, msa)
        self.__get_msa_position = get_msa_position

    def __record_to_position__(self, msa: MultipleSeqAlignment, record: TModelRecord) -> Optional[PositionEntry]:
        return self.__get_msa_position(msa, record)