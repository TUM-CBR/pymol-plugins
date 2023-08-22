from typing import Dict, Iterable, NamedTuple, Optional, Tuple

from .MeltingTemp import MeltingTemp

CODON_SIZE = 3

class PrimerResult(NamedTuple):
    left_primer : str
    tm_left : float
    right_primer : str
    tm_right : float
    inner_seq : str
    tm_all : float

    def __get_error(self, target_tm : float) -> float:
        return sum([
            abs(self.tm_left - target_tm),
            abs(self.tm_right - target_tm),
            abs(self.tm_all - target_tm)
        ])

    @staticmethod
    def choose_best(target_tm : float, opt1 : 'PrimerResult', opt2 : 'PrimerResult') -> 'PrimerResult':
        if opt1.__get_error(target_tm) < opt2.__get_error(target_tm):
            return opt1
        else:
            return opt2

class Operations:

    class DesignPrimers(NamedTuple):
        operations : 'Operations'
        sequence : str
        tm : float
        start : int
        count : int
        min_lenght : int
        max_length : int

        def design_primers(self) -> Iterable[Dict[str, PrimerResult]]:

            for i in range(self.start, self.count):
                yield self.design_primer_at(i)

        def design_primer_at(self, position : int) -> Dict[str, PrimerResult]:
            tm_calc = self.operations.tm_calc
            results = {}

            for (p_left, o_codon, p_right) in self.generate_primers_at(position):
                tm_left = tm_calc.oligo_tm(p_left)
                tm_right = tm_calc.oligo_tm(p_right)

                for aa,codons in CODONS.items():
                    codon = codons[0]
                    seq = p_left + o_codon + p_right
                    tm_all = tm_calc.oligo_tm_mis(seq, {len(seq): codon})
                    result = PrimerResult(
                        left_primer=p_left,
                        tm_left=tm_left,
                        right_primer=p_right,
                        tm_right=tm_right,
                        inner_seq=codon,
                        tm_all=tm_all
                    )

                    best_result = results.get(aa) or result
                    results[aa] = PrimerResult.choose_best(self.tm, best_result, result)

            return results

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

    def __init__(self, tm_calc : Optional[MeltingTemp]):
        self.__tm_calc = tm_calc or MeltingTemp()

    @property
    def tm_calc(self) -> MeltingTemp:
        return self.__tm_calc

    def design_primers(
        self,
        sequence : str,
        tm : float,
        start : int,
        count : int
    ) -> Iterable[Dict[str, PrimerResult]]:
        
        return Operations.DesignPrimers(
            sequence=sequence,
            tm=tm,
            start=start,
            count=count,
            min_lenght=6,
            max_length=36
        ).design_primers()


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
