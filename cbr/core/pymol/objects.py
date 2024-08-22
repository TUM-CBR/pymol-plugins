import pymol
from pymol import cmd, CmdException
from typing import Iterable, NamedTuple, Set, Tuple

class ModelAndChain(NamedTuple):
    model: str
    chain: str

    def iter_segs(self):
        segs: Set[str] = set()
        pymol.cmd.iterate(
            "model %s and chain %s" % (self.model, self.chain),
            'segs.add(segi)',
            space={'segs': segs}
        )

        return segs

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

def iter_segis() -> Iterable[Tuple[str, str, str]]:

    for (name, chain) in iter_chains():
        segs: Set[str] = set()
        pymol.cmd.iterate(
            "model %s and chain %s" % (name, chain),
            'segs.add(segi)',
            space={'segs': segs}
        )

        for seg in segs:
            yield (name, chain, seg)
