# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SchemaEnergyRunner.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SchemaEnergyRunner(object):
    def setupUi(self, SchemaEnergyRunner):
        SchemaEnergyRunner.setObjectName("SchemaEnergyRunner")
        SchemaEnergyRunner.resize(682, 518)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(SchemaEnergyRunner)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(SchemaEnergyRunner)
        font = QtGui.QFont()
        font.setPointSize(20)
        font.setBold(True)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(SchemaEnergyRunner)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.structuresCombo = QtWidgets.QComboBox(SchemaEnergyRunner)
        self.structuresCombo.setObjectName("structuresCombo")
        self.horizontalLayout_2.addWidget(self.structuresCombo)
        self.refreshStructuresButton = QtWidgets.QPushButton(SchemaEnergyRunner)
        self.refreshStructuresButton.setObjectName("refreshStructuresButton")
        self.horizontalLayout_2.addWidget(self.refreshStructuresButton)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.gridLayout.addLayout(self.horizontalLayout_2, 0, 0, 1, 1)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_4 = QtWidgets.QLabel(SchemaEnergyRunner)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_4.addWidget(self.label_4)
        self.structureSequenceCombo = QtWidgets.QComboBox(SchemaEnergyRunner)
        self.structureSequenceCombo.setObjectName("structureSequenceCombo")
        self.horizontalLayout_4.addWidget(self.structureSequenceCombo)
        self.gridLayout.addLayout(self.horizontalLayout_4, 1, 0, 1, 2)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.runSchemaEnergyButton = QtWidgets.QPushButton(SchemaEnergyRunner)
        self.runSchemaEnergyButton.setObjectName("runSchemaEnergyButton")
        self.horizontalLayout_5.addWidget(self.runSchemaEnergyButton)
        self.schemaProgress = QtWidgets.QProgressBar(SchemaEnergyRunner)
        self.schemaProgress.setProperty("value", 24)
        self.schemaProgress.setObjectName("schemaProgress")
        self.horizontalLayout_5.addWidget(self.schemaProgress)
        self.gridLayout.addLayout(self.horizontalLayout_5, 2, 2, 1, 1)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_5 = QtWidgets.QLabel(SchemaEnergyRunner)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_6.addWidget(self.label_5)
        self.shufflingPointsEdit = QtWidgets.QLineEdit(SchemaEnergyRunner)
        self.shufflingPointsEdit.setObjectName("shufflingPointsEdit")
        self.horizontalLayout_6.addWidget(self.shufflingPointsEdit)
        self.gridLayout.addLayout(self.horizontalLayout_6, 0, 2, 1, 1)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_6 = QtWidgets.QLabel(SchemaEnergyRunner)
        self.label_6.setObjectName("label_6")
        self.horizontalLayout_7.addWidget(self.label_6)
        self.energyScoringCombo = QtWidgets.QComboBox(SchemaEnergyRunner)
        self.energyScoringCombo.setObjectName("energyScoringCombo")
        self.horizontalLayout_7.addWidget(self.energyScoringCombo)
        self.gridLayout.addLayout(self.horizontalLayout_7, 1, 2, 1, 1)
        self.verticalLayout_2.addLayout(self.gridLayout)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_3 = QtWidgets.QLabel(SchemaEnergyRunner)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_3.addWidget(self.label_3)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem2)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.fastaTextEdit = QtWidgets.QPlainTextEdit(SchemaEnergyRunner)
        self.fastaTextEdit.setObjectName("fastaTextEdit")
        self.verticalLayout.addWidget(self.fastaTextEdit)
        self.verticalLayout_2.addLayout(self.verticalLayout)

        self.retranslateUi(SchemaEnergyRunner)
        QtCore.QMetaObject.connectSlotsByName(SchemaEnergyRunner)

    def retranslateUi(self, SchemaEnergyRunner):
        _translate = QtCore.QCoreApplication.translate
        SchemaEnergyRunner.setWindowTitle(_translate("SchemaEnergyRunner", "Form"))
        self.label.setText(_translate("SchemaEnergyRunner", "Schema Energy"))
        self.label_2.setText(_translate("SchemaEnergyRunner", "Select Structure"))
        self.refreshStructuresButton.setText(_translate("SchemaEnergyRunner", "Refresh"))
        self.label_4.setText(_translate("SchemaEnergyRunner", "Select Structure Sequence"))
        self.runSchemaEnergyButton.setText(_translate("SchemaEnergyRunner", "Run SCHEMA Energy"))
        self.label_5.setText(_translate("SchemaEnergyRunner", "Shuffling Points"))
        self.label_6.setText(_translate("SchemaEnergyRunner", "Energy Scoring"))
        self.label_3.setText(_translate("SchemaEnergyRunner", "Input Sequences (fasta)"))
