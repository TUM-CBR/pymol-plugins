from Bio.Align import MultipleSeqAlignment, SeqRecord

from typing import List

from ...clustal.msa import is_blank

def score_length(sequences : MultipleSeqAlignment) -> List[float]:

    def len_no_gaps(seq : SeqRecord) -> int:
        result = 0
        for resi in seq:
            if not is_blank(resi):
                result += 1

        return result

    return [len_no_gaps(seq) for seq in sequences]