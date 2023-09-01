from typing_extensions import Literal
from Bio.Seq import Seq
from typing import Dict, Iterable, List, NamedTuple, Optional, Tuple

from .MeltingTemp import MeltingTemp

CODON_SIZE = 3

PrimerOrganism = Literal['E_COLI', 'P_PASTORIS']

class PrimerResult(NamedTuple):
    left_primer : str
    tm_left : float
    right_primer : str
    tm_right : float
    inner_seq : str
    tm_all : Optional[float]
    tm_error : Optional[str]

    def __get_penalty(self, target_tm : float, primer_tm : float):
        return abs(target_tm - primer_tm)

    def __get_error(
        self,
        target_tm : float,
        w_target_tm : float,
        primers_tm : float,
        w_primers_tm : float,
        w_tm_delta : float
    ) -> float:

        return sum([
            w_primers_tm*self.__get_penalty(primers_tm, self.tm_left),
            w_primers_tm*self.__get_penalty(primers_tm, self.tm_right),
            w_tm_delta*abs(self.tm_left - self.tm_right),
            w_target_tm*(self.tm_all and abs(self.tm_all - target_tm) or 0)
        ])

    @staticmethod
    def choose_best(
        opt1 : 'PrimerResult',
        opt2 : 'PrimerResult',
        target_tm : float,
        w_target_tm : float,
        primers_tm : float,
        w_primers_tm : float,
        w_tm_delta : float
    ) -> 'PrimerResult':
        score_opt1 = opt1.__get_error(
            target_tm,
            w_target_tm,
            primers_tm,
            w_primers_tm,
            w_tm_delta
        )
        score_opt2 = opt2.__get_error(
            target_tm,
            w_target_tm,
            primers_tm,
            w_primers_tm,
            w_tm_delta
        )
        if score_opt2 < score_opt1:
            return opt2
        else:
            return opt1

    @property
    def c_left_primer(self):
        return str(Seq(self.left_primer).complement())

    @property
    def c_right_primer(self):
        return str(Seq(self.right_primer).complement())

    @property
    def c_inner(self):
        return str(Seq(self.inner_seq).complement())

DesignPrimersResult = Dict[str, List[PrimerResult]]
DesignPrimersResults = Iterable[DesignPrimersResult]

class Operations:

    class DesignPrimers(NamedTuple):
        operations : 'Operations'
        sequence : str
        start : int
        count : int
        min_lenght : int
        max_length : int
        codons : Dict[str, str]

        @property
        def __step(self):
            return CODON_SIZE

        def design_primers(self) -> DesignPrimersResults:

            for i in range(self.start, self.start + self.count, self.__step):
                yield self.design_primer_at(i)

        def design_primer_at(self, position : int) -> DesignPrimersResult:
            tm_calc = self.operations.tm_calc
            results : DesignPrimersResult = {}

            for aa,codon in self.codons.items():

                def get_primers_for_codon() -> Iterable[PrimerResult]:
                    for (p_left, o_codon, p_right) in self.generate_primers_at(position):
                        tm_left = tm_calc.oligo_tm(p_left)
                        tm_right = tm_calc.oligo_tm(p_right)

                        seq = p_left + o_codon + p_right
                        tm_error = None

                        try:
                            tm_all = tm_calc.oligo_tm_mis(seq, {len(p_left): codon})
                        except Exception as e:
                            tm_all = None
                            tm_error = str(e)
                        yield PrimerResult(
                            left_primer=p_left,
                            tm_left=tm_left,
                            right_primer=p_right,
                            tm_right=tm_right,
                            inner_seq=codon,
                            tm_all=tm_all,
                            tm_error=tm_error
                        )

                results[aa] = list(get_primers_for_codon())

            return results

        def __is_unique(self, primer_candidate: str):
            start_ix = 0

            while(start_ix < len(self.sequence)):

                next_ix = self.sequence.find(primer_candidate, start_ix)
                if next_ix >= 0 and start_ix > 0:
                    return False
                elif next_ix >= 0:
                    start_ix = next_ix + 1
                elif start_ix == 0:
                    raise Exception("The string %s does not occur in %s" % (primer_candidate, self.sequence))
                else:
                    return True

            return True

        def generate_primers_at(
            self,
            position : int,
        ) -> Iterable[Tuple[str, str, str]]:
            count = self.__step

            for i in range(self.min_lenght, self.max_length):
                for j in range(self.min_lenght, self.max_length):
                    l_start = position - i
                    r_start = position + count
                    r_end = position + count + j

                    if l_start > 0 and r_start < len(self.sequence):

                        l_seq = self.sequence[l_start:position]
                        r_seq = self.sequence[r_start:r_end]
                        if self.__is_unique(l_seq) and self.__is_unique(r_seq):
                            yield (
                                l_seq,
                                self.sequence[position:r_start],
                                r_seq
                            )

    def __init__(self, tm_calc : Optional[MeltingTemp] = None):
        self.__tm_calc = tm_calc or MeltingTemp()

    @property
    def tm_calc(self) -> MeltingTemp:
        return self.__tm_calc

    def design_primers(
        self,
        sequence : str,
        start : int,
        count : int,
        organism : PrimerOrganism
    ) -> DesignPrimersResults:

        return Operations.DesignPrimers(
            sequence=sequence,
            start=start,
            count=count,
            min_lenght=6,
            max_length=46,
            operations=self,
            codons=CODONS_MAP[organism]
        ).design_primers()

P_PASTORIS_CODONS = {
    "Ala": "GCT",
    "Cys": "TGT",
    "Asp": "GAC",
    "Glu": "GAG",
    "Phe": "TTC",
    "Gly": "GGT",
    "His": "CAC",
    "Ile": "ATT",
    "Lys": "AAG",
    "Leu": "TTG",
    "Met": "ATG",
    "Asn": "AAC",
    "Pro": "CCA",
    "Gln": "CAA",
    "Arg": "AGA",
    "Ser": "TCT",
    "Thr": "ACT",
    "Val": "GTT",
    "Trp": "TGG",
    "Tyr": "TAC"
}

E_COLI_CODONS = {
    "Ala": "GCG",
    "Cys": "TGC",
    "Asp": "GAC",
    "Glu": "GAA",
    "Phe": "TTC",
    "Gly": "GGT",
    "His": "CAC",
    "Ile": "ATC",
    "Lys": "AAA",
    "Leu": "CTG",
    "Met": "ATG",
    "Asn": "AAC",
    "Pro": "CCG",
    "Gln": "CAG",
    "Arg": "CGT",
    "Ser": "TCT",
    "Thr": "ACC",
    "Val": "GTT",
    "Trp": "TGG",
    "Tyr": "TAC"
}

CODONS_MAP : Dict[PrimerOrganism, Dict[str, str]] = {
    'E_COLI': E_COLI_CODONS,
    'P_PASTORIS': P_PASTORIS_CODONS
}