# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SchemaRunner.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SchemaRunner(object):
    def setupUi(self, SchemaRunner):
        SchemaRunner.setObjectName("SchemaRunner")
        SchemaRunner.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(SchemaRunner)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.label_2 = QtWidgets.QLabel(SchemaRunner)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.sequencesText = QtWidgets.QPlainTextEdit(SchemaRunner)
        self.sequencesText.setObjectName("sequencesText")
        self.gridLayout.addWidget(self.sequencesText, 2, 0, 1, 4)
        self.structuresCombo = QtWidgets.QComboBox(SchemaRunner)
        self.structuresCombo.setObjectName("structuresCombo")
        self.gridLayout.addWidget(self.structuresCombo, 0, 1, 1, 1)
        self.refreshButton = QtWidgets.QPushButton(SchemaRunner)
        self.refreshButton.setObjectName("refreshButton")
        self.gridLayout.addWidget(self.refreshButton, 0, 2, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 3, 1, 1)
        self.label = QtWidgets.QLabel(SchemaRunner)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.runSchemaButton = QtWidgets.QPushButton(SchemaRunner)
        self.runSchemaButton.setObjectName("runSchemaButton")
        self.gridLayout.addWidget(self.runSchemaButton, 4, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(SchemaRunner)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 3, 0, 1, 1)
        self.lineEdit = QtWidgets.QLineEdit(SchemaRunner)
        self.lineEdit.setObjectName("lineEdit")
        self.gridLayout.addWidget(self.lineEdit, 3, 1, 1, 3)
        self.verticalLayout.addLayout(self.gridLayout)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)

        self.retranslateUi(SchemaRunner)
        QtCore.QMetaObject.connectSlotsByName(SchemaRunner)

    def retranslateUi(self, SchemaRunner):
        _translate = QtCore.QCoreApplication.translate
        SchemaRunner.setWindowTitle(_translate("SchemaRunner", "Form"))
        self.label_2.setText(_translate("SchemaRunner", "Select Parents"))
        self.refreshButton.setText(_translate("SchemaRunner", "Refresh"))
        self.label.setText(_translate("SchemaRunner", "Select Structure:"))
        self.runSchemaButton.setText(_translate("SchemaRunner", "Run SCHEMA"))
        self.label_3.setText(_translate("SchemaRunner", "Corssover Counts"))
