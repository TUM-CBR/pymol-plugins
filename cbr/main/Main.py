from . import resources # pyright: ignore

from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QFont, QPixmap
from PyQt5.QtWidgets import QGridLayout, QHBoxLayout, QLabel, QSizePolicy, QSpacerItem, QWidget, QVBoxLayout
from os import path
from typing import Any, Callable, NamedTuple, Optional

from ..core.Context import Context
from .. import apbs as Apbs
from ..cascades import main as Cascades
from ..cavities import main as Cavitites
from ..chimeras import Applications as Chimeras
from ..coevolution import main as Coevolution
from ..dssp import main as Dssp
from ..frankenprot import main as FrankenProt
#from ..kinetics import main as Kinetics
from ..mpnn import main as Mpnn
from ..msa import Applications as Msa
from ..primer import Applications as Primer
from ..schema import Applications as Schema
from .AppIcon import AppIcon

class AppDefinition(NamedTuple):
    icon : str
    text : str
    init : Callable[[Context], Any]
    user_manual: Optional[QUrl] = None

APP_DEFINITIONS = [
    AppDefinition(
        icon = "schema.png",
        text = "SCHEMA RASPP",
        init = Schema.schema_raspp,
        user_manual = QUrl("https://github.com/TUM-CBR/pymol-plugins/wiki/SCHEMA-RASPP")
    ),
    AppDefinition(
        icon = "schema.png",
        text = "SCHEMA Energy",
        init = Schema.schema_energy,
        user_manual = QUrl("https://github.com/TUM-CBR/pymol-plugins/wiki/SCHEMA-Energy")
    ),
    AppDefinition(
        icon = "msa.png",
        text = "MSA Viewer",
        init = Msa.msa_viewer,
        user_manual = QUrl("https://github.com/TUM-CBR/pymol-plugins/wiki/MSA-Viewer")
    ),
    AppDefinition(
        icon = "msa.png",
        text = "MSA Cleaner",
        init = Msa.msa_cleaner
    ),
    AppDefinition(
        icon = "chimeras.jpg",
        text = "Chimeras Generator",
        init = Chimeras.chimeras_generator
    ),
    AppDefinition(
        icon = "primer.png",
        text = "Primers Designer",
        init = Primer.primer_generator
    ),
    AppDefinition(
        icon = "coevolution.png",
        text = "Coevolution Analysis",
        init = Coevolution.coevolution_runner
    ),
    AppDefinition(
        icon = "apbs.png",
        text = "APBS Electrostatics",
        init = Apbs.run_apbs_electrostatics,
        user_manual = QUrl("https://pymolwiki.org/index.php/APBS")
    ),
    AppDefinition(
        icon = "cascades.png",
        text = "Cascade BLAST",
        init = Cascades.cascades_main,
        user_manual=QUrl("https://github.com/TUM-CBR/pymol-plugins/wiki/Cascade-BLAST")
    ),
    AppDefinition(
        icon = "dssp-logo.svg",
        text = "DSSP",
        init = Dssp.dssp_runner
    ),
    AppDefinition(
        icon = "proteinMPNN.png",
        text = "ProteinMPNN",
        init = Mpnn.mpnn_main,
        user_manual = QUrl('https://github.com/TUM-CBR/pymol-plugins/wiki/ProteinMPNN')
    ),
    AppDefinition(
        icon = "fasta.png",
        text = "Fasta Viewer",
        init = Msa.fasta_viewer
    ),
    AppDefinition(
        icon = "frankenprot.png",
        text = "FrankenProt",
        init = FrankenProt.frankenprot
    ),
    AppDefinition(
         icon = "protodentist.webp",
         text = "ProtoDentist",
         init = Cavitites.main
    )
    #AppDefinition(
    #    icon = "kinetics.png",
    #    text = "Enzyme Kinetics",
    #    init = Kinetics.enzyme_kinetics
    #)
]

class Main(QWidget):

    def __init__(self, context : Context, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

        self.__context = context
        self.__init_ui()

    def __get_image_path(self, name : str) -> str:
        return path.join(
            path.dirname(__file__),
            "resources",
            name
        )

    def __init_ui(self):

        main_layout = QVBoxLayout()

        title_layout = QHBoxLayout()

        pixmap = QPixmap(self.__get_image_path("tumcs.png"))
        pixmap = pixmap.scaledToWidth(160, Qt.TransformationMode.SmoothTransformation)
        image_label = QLabel()
        image_label.setPixmap(pixmap)
        title_layout.addWidget(image_label)

        label = QLabel()
        font = QFont()
        font.setPointSize(28)
        font.setBold(True)
        label.setFont(font)
        label.setText("CBR Bioinformatics Tools")

        title_layout.addWidget(label)
        title_layout.addItem(
            QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        )

        main_layout.addLayout(title_layout)

        apps_layout = QGridLayout()

        for i,app_definition in enumerate(APP_DEFINITIONS):
             self.__add_app(i, app_definition, apps_layout)

        main_layout.addLayout(apps_layout)
        main_layout.addItem(
            QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        )

        self.setLayout(main_layout)


    def __add_app(self, i: int, app_definition : AppDefinition, apps_layout : QGridLayout):
            columns = 4
            row = int(i / columns)
            col = i % columns
            icon = self.__get_image_path(app_definition.icon)

            def __init_app__():
                app_definition.init(
                    self.__context.create_child_context(
                        new_user_manual=app_definition.user_manual # type: ignore
                    )
                )

            app_icon = AppIcon(app_definition.text, icon, __init_app__)
            apps_layout.addWidget(app_icon, row, col)
