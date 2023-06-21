import pymol.cmd as cmd
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QFileDialog, QWidget
from typing import Dict, Tuple

from ...core.Context import Context
from ...core import visual
from ...core.pymol import structure
from ...clustal import msa
from ...clustal import Clustal
from .Ui_MsaViewer import Ui_MsaViewer

COLOR_MAX = 999
COLOR_MIN = 600

class MsaViewer(QWidget):

    def __init__(self, context : Context, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.__ui = Ui_MsaViewer()
        self.__ui.setupUi(self)
        visual.as_structure_selector(self.__ui.structuresCombo, self.__ui.structuresRefreshButton)
        self.__set_sequences({})
        self.__clustal = Clustal.get_clustal_from_context(context)

    def __set_sequences(self, sequences : Dict[str, str]):
        self.__sequences = dict(sequences)
        self.__ui.sequenceCombo.clear()
        self.__ui.sequenceCombo.addItems(self.__sequences.keys())

    @pyqtSlot()
    def on_selectMsaButton_clicked(self):
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly  # Optionally set file dialog options

        # Show the file dialog and get the selected file
        file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "FASTA alignment Files (*.*)", options=options)

        if not file_path:
            return

        print("Selected file:", file_path)
        self.__set_sequences(msa.parse_alignments(file_path))

    @property
    def __selected_structure(self) -> Tuple[str, str]:
        return self.__ui.structuresCombo.currentData()

    def __get_structure_sequence(self):
        (structure_name, chain) = self.__selected_structure
        return structure.get_pdb_sequence(structure_name, chain)

    def __get_structure_positions(self):
        sequence_name = self.__ui.sequenceCombo.currentText()
        sequence = msa.clean_msa_blanks(self.__sequences[sequence_name])
        structure_sequence = self.__get_structure_sequence()
        result = self.__clustal.run_msa_items(
            [ ("structure", structure_sequence)
            , (sequence_name, sequence)
            ]
        )

        return list(msa.get_relative_positions(self.__sequences, result))

    def __get_structure_conservation(self):
        sequences = list(self.__sequences.values())
        positions = self.__get_structure_positions()
        result = {}

        for i in range(0, len(sequences[0])):
            ix = positions[i]
            conserved = {}
            for sequence in sequences:
                aa = sequence[i].lower()
                if msa.is_blank(aa):
                    continue

                if aa in conserved:
                    conserved[aa] += 1
                else:
                    conserved[aa] = 1

            combined = sum(2 ** x for x in conserved.values())
            #combined = max(conserved.values())
            if ix in result:
                result[ix] = max(combined, result[ix])
            else:
                result[ix] = combined

        return result

    def __color_by_conservation(self):
        (selected_structure, chain) = self.__selected_structure
        offset = structure.get_structure_offset(selected_structure, chain)
        conservation = self.__get_structure_conservation()
        factor = (COLOR_MAX - COLOR_MIN) / max(conservation.values())
        structure_query = structure.get_structure_query(selected_structure, chain) 

        cmd.color("999", structure_query)

        for (i,score) in conservation.items():
            i_score = int(score * factor)
            ix = i + offset
            cmd.color(str(999 - i_score), "%s and resi %i-%i" % (structure_query, ix, ix))

    @pyqtSlot()
    def on_colorConservedButton_clicked(self):
        self.__color_by_conservation()
