from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from PyQt5.QtCore import QObject
from PyQt5.QtGui import QColor
from typing import Callable, Dict, Iterable, List, NamedTuple, Optional, Sequence

from ...clustal.Clustal import Clustal
from ...core.extensions import biopython as BioExtra
from ...core.pymol.structure import StructureSelection
from ...support.display.sequence import RESIDUE_COLORS

from .AbstractCompositeTableModel import *
from .StructuresAlignmentMapper import StructureAlignmentEntry, StructureResidue, StructuresAlignmentMapper

class SetStructureArg(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str

class PositionView(NamedTuple):
    positions: Sequence[int]
    structure_color: Optional[QColor]
    residues: Sequence[str]

class PositionsWithAttributes(NamedTuple):
    positions: Sequence[int]
    structure_color: Optional[QColor] = None

    def as_view(self, msa_mapping: Dict[int, StructureResidue]) -> Optional['PositionView']:

        mapped_positions = [msa_mapping[p] for p in self.positions if p in msa_mapping]

        if len(mapped_positions) == 0:
            return None
        else:
            return PositionView(
                positions = [p.resv for p in mapped_positions],
                structure_color = self.structure_color, # type: ignore
                residues = [p.resv_name for p in mapped_positions]
            )

class PositionEntry(NamedTuple):
    positions_with_attributes: Sequence[PositionsWithAttributes]

    @property
    def positions(self) -> Sequence[int]:
        return [
            position
            for entry in self.positions_with_attributes
            for position in entry.positions
        ]
    
    def map_positions(self, msa_mapping: Dict[int, StructureResidue]) -> Iterable[PositionView]:

        return (
            mapped
            for p in self.positions_with_attributes
            for mapped in [p.as_view(msa_mapping)]
            if mapped is not None
        )
    
    def view_msa(
        self,
        consensus: Seq
    ) -> List[ViewRecordAttributes]:
        
        positions = [
            pos
            for attributes in self.positions_with_attributes
            for pos in attributes.positions
        ]

        names = [
            consensus[p]
            for p in positions
        ]

        if len(names) == 1:
            background_color = RESIDUE_COLORS[names[0]]
        else:
            background_color = None

        return [
            ViewRecordAttributes(
                display=",".join(str(p) for p in positions)
            ),
            ViewRecordAttributes(
                display=",".join(n for n in names),
                background_color=background_color
            )
        ]
    
    def view_structure(
        self,
        msa_mapping: Dict[int, StructureResidue]
    ) -> List[ViewRecordAttributes]:

        records = list(self.map_positions(msa_mapping))

        positions = [
            p
            for mapped in records
            for p in mapped.positions
        ]

        names = [
            p
            for mapped in records
            for p in mapped.residues
        ]

        assert len(positions) == len(names)

        return [
            ViewRecordAttributes(
                display=",".join(str(p) for p in positions)
            ),
            ViewRecordAttributes(
                display=",".join(n for n in names)
            )
        ]
    
    def as_pymol_attribute(
        self,
        structure: StructureSelection,
        msa_mapping: Dict[int, StructureResidue]
    ) -> List[PymolRecordAttributes]:

        selections: List[StructureSelection] = []
        colors: List[Optional[QColor]] = []

        for mapped in self.map_positions(msa_mapping):
            selections.append(structure.scoped(mapped.positions))
            colors.append(mapped.structure_color) # type: ignore

        assert len(selections) == len(colors), "Algorithm is inconsistent"

        if len(selections) == 0:
            attributes = [PymolRecordAttributes.default()] * AbstractStructurePositionView.attributes_per_record()
        else:
            attributes = [
                PymolRecordAttributes(
                    selections=selections,
                    colors=colors
                ),
                PymolRecordAttributes.default()
            ]

        assert len(attributes) == AbstractStructurePositionView.attributes_per_record()

        return attributes

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
    __msa_to_resv: Sequence[Dict[int, StructureResidue]]

    K_POSITION = "Position"
    K_RES_NAME = "Residue"
    K_RECORD_HEADERS = [K_POSITION, K_RES_NAME]

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

    @classmethod
    def attributes_per_record(cls):
        return len(cls.K_RECORD_HEADERS)

    def __headers(self) -> Sequence[ViewHeaderSpec]:
        structure_headers = [
            f"{structure_header} ({structure.structure.show()})"
            for structure in self.__structures
            for structure_header in self.K_RECORD_HEADERS
        ]

        msa_headers = [
            f"{self.K_POSITION} (msa)",
            f"{self.K_RES_NAME} (msa consensus)"
        ]

        return [
            ViewHeaderSpec(name)
            for name in msa_headers + structure_headers
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

        self.content_changed.emit(self)

    def __record_to_position__(self, msa: MultipleSeqAlignment, record: TModelRecord) -> Optional[PositionEntry]:
        raise NotImplementedError("The function __record_to_position__ must be implemented.")
    
    def __get_record_attributes(
        self,
        mapper: StructuresAlignmentMapper,
        consensus: Seq,
        record: TModelRecord
    ) -> Sequence[ViewRecordAttributes]:

        position = self.__record_to_position__(mapper.msa, record)

        if position is None:
            position_view = self.DEFAULT_QT_ATTRIBUTES
            structures_view = [self.DEFAULT_QT_ATTRIBUTES for _ in self.__msa_to_resv]
            return [position_view] + structures_view

        position_view = position.view_msa(consensus)

        structures_view = [
            attr
            for mapping in self.__msa_to_resv
            for attr in position.view_structure(mapping)
        ]

        return position_view + structures_view
    
    def __get_pymol_attributes(self, mapper: StructuresAlignmentMapper, record: TModelRecord) -> Sequence[PymolRecordAttributes]:

        position = self.__record_to_position__(mapper.msa, record)

        default_record_attributes = [self.DEFAULT_PYMOL_ATTRIBUTES] * self.attributes_per_record()

        if position is None:
            return default_record_attributes \
                + [
                    attr
                    for _ in self.__msa_to_resv
                    for attr in default_record_attributes
                ]

        return default_record_attributes + \
            [
                attr
                for structure, mapping in zip(self.__structures, self.__msa_to_resv)
                for attr in position.as_pymol_attribute(structure.structure, mapping)
            ]

    def attributes(self, records: Sequence[TModelRecord]) -> ViewRecords:
        
        mapper = self.__mapper

        if mapper is None:
            return super().attributes(records)
        
        qt_attributes: Sequence[Sequence[ViewRecordAttributes]] = []
        pymol_attributes: Sequence[Sequence[PymolRecordAttributes]] = []
        consensus = BioExtra.consensus(mapper.msa)

        for record in records:
            qt_attributes.append(self.__get_record_attributes(mapper, consensus, record))
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