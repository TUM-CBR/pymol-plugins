# CBR Bioinformatics Tools
> Enjoyable Bioinformatics

CBR Bioinformatics Tools is a suite of software applications which allow researchers to use various structural and homology related algorithms and models with an intuitive graphical user interface.

The goal of this project is to make state of the art bioinformatics algorithms and models accessible to everyone, in a hassle free fashion requiring no configuration or programing knowledge. The following steps are taken to meet this goal:

* The suite is distributed as a [PyMol](https://pymol.org/2/) plug-in to allow easy integration with molecular visualization software.
* Instead of distributing standalone scripts, the suite si packaged along wiht all necessary dependencies required for proper operation.
* Some components are distributed separately (to keep the size in check). In such cases, a graphical installer automatically installs the features upon first use of the respective tool.
* Our packages are tested against the official PyMol release on the Windows operating system to ensure compatibility.

We hope you enjoy using our Bioinformatics suite. Please let us know if there is anything we can make to improve your experience by [creating an issue](https://github.com/TUM-CBR/pymol-plugins/issues).

## Table of Contents
- [Installation and Useage](#installation-and-useage)
- [Features](#features)
- [License](#license)
- [Credits](#credits)

## Installation and Useage

For information about how to install and use the suite, please refer to the [The CBR Bioinformatics Tools Wiki](https://github.com/TUM-CBR/pymol-plugins/wiki).

## Features

Information about the applications and features of this suite can be found in the [The CBR Bioinformatics Tools Wiki](https://github.com/TUM-CBR/pymol-plugins/wiki). We will briefly highlight some of the features below:

* [MSA Viewer](https://github.com/TUM-CBR/pymol-plugins/wiki/MSA-Viewer) Is a tool that allows coloring protein structures using statistical information from sequence alignment. Think of it as a "Color By MSA" feature.
* [Prot 'o Dentist](https://github.com/TUM-CBR/pymol-plugins/wiki/Prot-o'-Dentist) is a tool for identifying cavities in a protein and easily interacting with the residues surrounding the cavities.
* [ProteinMPNN](https://github.com/TUM-CBR/pymol-plugins/wiki/ProteinMPNN) is a graphical user frontend for the [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) AI model.
* [FrankenProt](https://github.com/TUM-CBR/pymol-plugins/wiki/FrankenProt) Is a tool that allows modifying protein backbones. It is meant to be used alongside [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) to generate sequences that will fold slightly differently than the template backbone.
* [Cascade BLAST](https://github.com/TUM-CBR/pymol-plugins/wiki/Cascade-BLAST) is a tool that uses BLAST to find organisms having similar enzymatic pathwas as a given input.

## License

CBR Bioinformatics Tools is released under the "GNU General Public License Version 2". This license allows you to freely use and re-distribtue the software (even for commercial and for-profit reasons) with the condition that any modifications or improvements are made available to the general public under the same license.

Furthermore, these tools depend on other programs to operate correctly. In particular:

* [PyMOL](https://pymol.org/2/)
* [QT](https://www.qt.io/product)
* [Clustal Omega](http://www.clustal.org/clustal2/)
* [APBS](https://apbs.readthedocs.io/en/latest/)
* [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)

You may need additional licenses from these vendors to use and/or distribute this suite in a commercial fashion. Please refer to their sites for more information.

## Credits

This suite is develooped by the [Lehrstuhl für Chemie Biogener Rohstoffe](https://cbr.cs.tum.de/) with funding from the [Technical University of Munich](https://www.tum.de/) and the [State of Bavaria](https://www.bayern.de/).

Special acknowledgements for the team directly involved in desiging, developing and testing this suite:
* [Prof. Dr. Volker Sieber](https://cbr.cs.tum.de/team/volker-sieber/)
* [Dr. Amelie Skopp](https://cbr.cs.tum.de/team/amelie-skopp-dr/)
* [Dr. Enrico Hupfeld](https://cbr.cs.tum.de/team/enrico-hupfeld/)
* [MsC. Ernesto Rodriguez](https://github.com/netogallo)