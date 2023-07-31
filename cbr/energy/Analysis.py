import csv

from ..core.sequence import aa_from_letter
from .Foundations import GromacsRun
from .ProThermDB import ProThermEntry
from typing import Iterable, List, NamedTuple, Optional, TextIO, Tuple

class PotentialCsvRow(NamedTuple):
    structure : str
    position : Optional[int]
    mutation : Optional[str]
    energy : float
    ddG : Optional[float]
    pH : Optional[str]

def write_potential_csv(
    items : Iterable[GromacsRun],
    out_stream : TextIO,
    pro_therm_entries : Optional[Iterable[ProThermEntry]] = None
    ) -> None:
    
    pro_therm  : dict[Tuple[str, Optional[int], str], List[ProThermEntry]] = {}
    for entry in pro_therm_entries or []:

        mutation = entry.mutation
        if mutation is None:
            continue

        k = (entry.pdb_code.lower(), mutation.position, aa_from_letter(mutation.new_residue).lower())
        if k not in pro_therm:
            pro_therm[k] = []

        pro_therm[k].append(entry)

    csv_writer = csv.writer(out_stream)

    for item in items:
        structure = item.metadata.name.split("_")[0]
        mutation = next(item.metadata.mutations.__iter__(), None)

        if mutation:
            position = int(mutation.selection.replace("resi", ""))
            new_residue = item.metadata.mutations[0].mutation
        else:
            position = None
            new_residue = None

        try:
            energy_entry = max(
                item.emin_log or [],
                key = lambda i: i.step or -1
            )
        except StopIteration:
            print("Spurious entry %s" % item.metadata.name)
            continue

        pro_therm_items = pro_therm.get((structure.lower(), position, (new_residue or "").lower()))

        row = PotentialCsvRow(
            structure = structure,
            position = position,
            mutation = new_residue,
            energy = energy_entry.potential,
            ddG = None,
            pH = None
        )

        if pro_therm_items:
            for item in pro_therm_items:
                if item.ddG:
                    csv_writer.writerow(row._replace(ddG = item.ddG, pH = item.pH))
        else:
            csv_writer.writerow(row)

