from PyQt5.QtWidgets import QComboBox
from os import path
from typing import Optional

from ..raspp import contacts

class EnergySelector(object):

    interactions = {
        'SCHEMA classic' : None,
        'Simplified Pysics' : [
            contacts.electrostatic_interactions,
            contacts.van_der_waals,
            contacts.h_bond_interactions
        ]
    }

    def __init__(self, energyScoringCombo : QComboBox):
        self.__energyScoringCombo = energyScoringCombo
        for (i, (name, data)) in enumerate(EnergySelector.interactions.items()):
            energyScoringCombo.insertItem(
                i,
                name,
                userData=data
            )

    def __schema_interactions_file(self, base_path : str) -> str:
        return path.join(
            base_path,
            "schema_interactions.json"
        )

    def write_interactions(self, base_path : str) -> Optional[str]:
        scoring = self.__energyScoringCombo.currentData()

        if scoring is None:
            return None
        else:
            f_name = self.__schema_interactions_file(base_path)
            with open(f_name, 'w') as f_interactions:
                contacts.write_interactions(scoring, f_interactions)
            return f_name