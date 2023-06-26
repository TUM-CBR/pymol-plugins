from PyQt5.QtWidgets import QTableWidgetItem, QWidget
from pymol import cmd

from ...core.Context import Context
from ..support import parsers
from .Ui_SchemaEnergyViewer import Ui_SchemaEnergyViewer


class SchemaEnergyViewer(QWidget):

    def __init__(
        self,
        context : Context,
        structure_name : str,
        chain_name : str,
        raw_results : str,
        raw_contacts : str,
        *args,
        **kwargs):
        super(SchemaEnergyViewer, self).__init__(*args, **kwargs)
        self.__ui = Ui_SchemaEnergyViewer()
        self.__ui.setupUi(self)
        self.__results = parsers.parse_schema_energy(raw_results)
        self.__contacts = parsers.parse_schema_contacts(raw_contacts)
        self.__structure_name = structure_name
        self.__chain_name = chain_name
        self.__render_everything()

    def __render_everything(self):

        results_table = self.__ui.resultsTable
        results_table.reset()
        entries = list(self.__results.entries)
        results_table.setRowCount(len(entries) + 1)

        results_table.setItem(0, 0, QTableWidgetItem('<average>'))
        results_table.setItem(0, 1, QTableWidgetItem(str(self.__results.average_energy)))
        results_table.setItem(0, 2, QTableWidgetItem(str(self.__results.average_mutations)))

        for (i,result) in enumerate(entries):
            # First row is the average
            i = i + 1
            results_table.setItem(i, 0, QTableWidgetItem(str(result.chimera)))
            results_table.setItem(i, 1, QTableWidgetItem(str(result.energy)))
            results_table.setItem(i, 2, QTableWidgetItem(str(result.mutations)))

        contacts_table = self.__ui.contactsTable
        contacts_table.reset()
        contacts = list(self.__contacts.contacts)
        contacts_table.setRowCount(len(contacts))
        residues = {}
        cmd.iterate(
            'model %s & chain %s' % (self.__structure_name, self.__chain_name),
            'residues[to_int(resi)] = resn',
            space = { 'residues' : residues, 'to_int': int }
        )

        for (i, contact) in enumerate(contacts):
            contacts_table.setItem(i, 0, QTableWidgetItem(str(contact.i)))
            contacts_table.setItem(i, 1, QTableWidgetItem(str(contact.j)))
            
            pdb_i = contact.pdb_i
            pdb_i_name = residues.get(pdb_i) or "?"
            contacts_table.setItem(
                i,
                2,
                QTableWidgetItem("%i (%s)" % (pdb_i, pdb_i_name))
            )

            pdb_j = contact.pdb_j
            pdb_j_name = residues.get(pdb_j) or "?"
            contacts_table.setItem(
                i,
                3,
                QTableWidgetItem("%i (%s)" % (pdb_j, pdb_j_name))
            )


