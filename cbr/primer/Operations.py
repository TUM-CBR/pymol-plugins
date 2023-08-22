from typing import Iterable, NamedTuple, Optional, Sequence, Tuple

from .Primer3 import Primer3

CODON_SIZE = 3

class PrimerResult(NamedTuple):
    left_primer : str
    tm_left : float
    right_primer : str
    tm_right : float
    inner_seq : str
    tm_all : float

    @staticmethod
    def __get_error(target_tm : float, result : 'PrimerResult') -> float:
        return sum([
            abs(result.tm_left - target_tm),
            abs(result.tm_right - target_tm),
            abs(result.tm_all - target_tm)
        ])

    def is_worse_than(self, target_tm : float, other : 'PrimerResult'):
        return self.__get_error(target_tm, self) > self.__get_error(target_tm, other)

class Operations:

    class DesignPrimers(NamedTuple):
        operations : 'Operations'
        sequence : str
        tm : float
        start : int
        count : int
        min_lenght : int
        max_length : int

        def design_primers(self):
            pass

        def design_primer_at(self, position : int):
            primer3 = self.operations.primer3
            results = {}

            for (p_left, _, p_right) in self.generate_primers_at(position):
                tm_left = primer3.oligo_tm(p_left)
                tm_right = primer3.oligo_tm(p_right)

                for aa,codons in CODONS.items():
                    codon = codons[0]
                    tm_all = primer3.oligo_tm(p_left + codon + p_right)
                    result = PrimerResult(
                        left_primer=p_left,
                        tm_left=tm_left,
                        right_primer=p_right,
                        tm_right=tm_right,
                        inner_seq=codon,
                        tm_all=tm_all
                    )


            pass

        def generate_primers_at(
            self,
            position : int,
        ) -> Iterable[Tuple[str, str, str]]:
            count = CODON_SIZE

            for i in range(self.min_lenght, self.max_length):
                for j in range(self.min_lenght, self.max_length):
                    l_start = position - i
                    r_start = position + count + j

                    if l_start > 0 and r_start < len(self.sequence):
                        yield (
                            self.sequence[l_start:position],
                            self.sequence[position:count],
                            self.sequence[count:r_start]
                        )

    def __init__(self, primer3 : Optional[Primer3]):
        self.__primer3 = primer3 or Primer3()

    @property
    def primer3(self) -> Primer3:
        return self.__primer3

    def design_primers(
        self,
        sequence : str,
        tm : float,
        start : int,
        count : int
    ):

        pass

    def __design_primers_at_position(
        self,
        args : DesignPrimersArgs,
        position : int
    ):
        

        pass

CODONS = {
  "Alanine": ["GCT", "GCC", "GCA", "GCG"],
  "Arginine": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
  "Asparagine": ["AAT", "AAC"],
  "Aspartic Acid": ["GAT", "GAC"],
  "Cysteine": ["TGT", "TGC"],
  "Glutamine": ["CAA", "CAG"],
  "Glutamic Acid": ["GAA", "GAG"],
  "Glycine": ["GGT", "GGC", "GGA", "GGG"],
  "Histidine": ["CAT", "CAC"],
  "Isoleucine": ["ATT", "ATC", "ATA"],
  "Leucine": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
  "Lysine": ["AAA", "AAG"],
  "Methionine": ["ATG"],
  "Phenylalanine": ["TTT", "TTC"],
  "Proline": ["CCT", "CCC", "CCA", "CCG"],
  "Serine": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
  "Threonine": ["ACT", "ACC", "ACA", "ACG"],
  "Tryptophan": ["TGG"],
  "Tyrosine": ["TAT", "TAC"],
  "Valine": ["GTT", "GTC", "GTA", "GTG"],
  "Stop Codon": ["TAA", "TAG", "TGA"]
}
