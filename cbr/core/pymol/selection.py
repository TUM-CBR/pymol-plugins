from typing import cast, List, Optional


class PymolSelection(object):

    def __init__(self, structrue : str, chain : Optional[str] = None):
        self.__structure = structrue
        self.__chain = chain

    def query(self) -> str:
        
        items = \
            list(
                filter(
                    lambda x: x,
                    [ "model %s" % self.__structure
                    , self.__chain and "chain %s" % self.__chain
                    ]
                )
            )
        return " & ".join(cast(List[str], items))

def get_selection_for_model(model : str, chain : Optional[str] = None):
    return PymolSelection(model, chain)