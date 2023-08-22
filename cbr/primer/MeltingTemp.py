from io import StringIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as BioTm
from typing import Dict, Optional

class MeltingTemp:

    def __exception(self, text : str) -> Exception:
        return Exception(text)

    def __oligo_tm_nn(self, seq : str, c_seq : Optional[str] = None) -> float:
        return BioTm.Tm_NN(seq, c_seq=c_seq)
    

    def oligo_tm(self, seq : str) -> float:
        return self.__oligo_tm_nn(seq)
    
    def oligo_tm_mis(self, seq : str, mismatches : Dict[int, str]) -> float:

        with StringIO(seq) as m_seq:

            for i,s in mismatches.items():
                m_seq.seek(i)
                m_seq.write(s)

            m_seq.seek(0)
            c_seq = str(Seq(m_seq.read()).complement())

        return self.__oligo_tm_nn(seq, c_seq=c_seq)