from PyQt5.QtCore import QUrl, pyqtSlot
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QTableWidgetItem, QWidget
from pymol import cmd
from tempfile import TemporaryDirectory

from ...core.Context import Context
from ..support import parsers
from .Ui_SchemaEnergyViewer import Ui_SchemaEnergyViewer


class SchemaEnergyViewer(QWidget):

    CONTACTS_HEADERS = [
        "i",
        "j",
        "pdb i",
        "pdb j",
    ]

    def __init__(
        self,
        context : Context,
        structure_name : str,
        chain_name : str,
        raw_results : str,
        raw_contacts : str,
        results_folder : TemporaryDirectory,
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
        self.__results_folder = results_folder

    def __del__(self):
        self.__results_folder.cleanup()

    @pyqtSlot()
    def on_openResultsFolder_clicked(self):
        QDesktopServices.openUrl(QUrl.fromLocalFile(self.__results_folder.name))

    @pyqtSlot()
    def on_contactsTable_itemSelectionChanged(self):

        selected = self.__ui.contactsTable.selectedIndexes()
        contacts = self.__contacts
        residues = " | ".join(
            "resi %i | resi %i" % (contacts.contacts[sel.row()].pdb_i, contacts.contacts[sel.row()].pdb_j)
            for sel in selected
        )

        if(len(residues) > 0):
            cmd.select(
                'SCHEMA_CONTACTS',
                "model %s & chain %s & (%s)" % (self.__structure_name, self.__chain_name, residues)
            )

    def __render_everything(self):

        results_table = self.__ui.resultsTable
        results_table.reset()
        entries = list(self.__results.entries)
        results_table.setRowCount(len(entries) + 1)

        interactions_meta = self.__contacts.interactions

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
        contacts = self.__contacts
        contacts_table.setRowCount(len(contacts.contacts))
        contacts_table.setColumnCount(4 + len(interactions_meta))
        contacts_table.setHorizontalHeaderLabels(
            SchemaEnergyViewer.CONTACTS_HEADERS + \
            [interaction.name for interaction in interactions_meta]
        )
        residues = {}
        cmd.iterate(
            'model %s & chain %s' % (self.__structure_name, self.__chain_name),
            'residues[to_int(resi)] = resn',
            space = { 'residues' : residues, 'to_int': int }
        )

        for (i, contact) in enumerate(contacts.contacts):
            contacts_table.setItem(i, 0, QTableWidgetItem(str(contact.seq_i)))
            contacts_table.setItem(i, 1, QTableWidgetItem(str(contact.seq_j)))
            
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

            interactions_summary = dict((meta.name, 0.0) for meta in interactions_meta)
            for interaction in contact.interactions:
                interactions_summary[interaction.name] += interaction.strength

            for (j, meta) in enumerate(interactions_meta):
                value = round(interactions_summary[meta.name], ndigits=2)
                contacts_table.setItem(i, j + 4, QTableWidgetItem(str(value)))


