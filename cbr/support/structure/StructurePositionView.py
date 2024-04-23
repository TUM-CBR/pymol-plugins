from Bio.Align import MultipleSeqAlignment
from PyQt5.QtCore import QObject
from PyQt5.QtGui import QColor
from typing import Dict, NamedTuple, Optional, Sequence

from cbr.support.structure.AbstractCompositeTableModel import PymolRecordAttributes, ViewHeaderSpec, ViewRecordAttributes

from ...clustal.Clustal import Clustal
from ...core.pymol.structure import StructureSelection

from .AbstractCompositeTableModel import AbstractRecordView, TModelRecord, ViewHeaderSpec
from .StructuresAlignmentMapper import StructureAlignmentEntry, StructuresAlignmentMapper

class SetStructureArg(NamedTuple):
    structure: StructureSelection
    reference_sequence_id: str

class PositionEntry(NamedTuple):
    position: int
    structure_color: Optional[QColor] = None

class AbstractMsaRecordView(AbstractRecordView[TModelRecord]):
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
        parent: Optional[QObject] = None
    ) -> None:
        super().__init__(parent)

        self.__clustal = clustal
        self.__mapper = None
        self.__structures = []
        self.__msa_to_resv = []

    def headers(self) -> Sequence[ViewHeaderSpec]:
        structure_headers = [
            f"Position ({structure.structure.show()})"
            for structure in self.__structures
        ]

        return [
            ViewHeaderSpec(name)
            for name in [self.K_POSITION] + structure_headers
        ]

    def __set_structures(self, structures: Sequence[SetStructureArg]):

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

        self.headers_changed.emit()

    def set_alignment(self, msa: MultipleSeqAlignment):

        self.__mapper = StructuresAlignmentMapper(
            msa,
            self.__clustal
        )

        self.__set_structures([
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

        position_view = ViewRecordAttributes(display=position.position)

        structures_view = [
            self.DEFAULT_QT_ATTRIBUTES if resv is None else ViewRecordAttributes(display=resv)
            for mapping in self.__msa_to_resv
            for resv in [mapping.get(position.position)]
        ]

        return [position_view] + structures_view

    def qt_attributes(self, records: Sequence[TModelRecord]) -> Sequence[Sequence[ViewRecordAttributes]]:

        mapper = self.__mapper

        if mapper is None:
            return super().qt_attributes(records)

        return [
            self.__get_record_attributes(mapper, record)
            for record in records
        ]
    
    def __get_pymol_attributes(self, mapper: StructuresAlignmentMapper, record: TModelRecord) -> Sequence[PymolRecordAttributes]:

        position = self.__record_to_position__(mapper.msa, record)

        if position is None:
            return [self.DEFAULT_PYMOL_ATTRIBUTES] \
                + [self.DEFAULT_PYMOL_ATTRIBUTES for _ in self.__msa_to_resv]

        return [self.DEFAULT_PYMOL_ATTRIBUTES] + \
            [
                self.DEFAULT_PYMOL_ATTRIBUTES \
                    if resv is None \
                    else PymolRecordAttributes(
                        selections=[structure.structure.scoped([resv])],
                        colors=[position.structure_color] # type: ignore
                    )
                for structure, mapping in zip(self.__structures, self.__msa_to_resv)
                for resv in [mapping.get(position.position)]
            ]

    def pymol_attributes(self, records: Sequence[TModelRecord]) -> Sequence[Sequence[PymolRecordAttributes]]:
        
        mapper = self.__mapper

        if mapper is None:
            return super().pymol_attributes(records)
        
        return [
            self.__get_pymol_attributes(mapper, record)
            for record in records
        ]