from pymol import cmd, CmdException
from typing import NamedTuple

class ModelAndChain(NamedTuple):
    model: str
    chain: str

def iter_chains():

    for model in cmd.get_names():
        try:
            for chain in cmd.get_chains(model):
                yield ModelAndChain(
                    model,
                    chain # type: ignore
                )
        except CmdException:
            # Some public objects throw errors when
            # iterating their chains
            pass