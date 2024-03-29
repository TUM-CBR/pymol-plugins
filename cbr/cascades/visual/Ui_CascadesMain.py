# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/cascades/visual/CascadesMain.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_CascadesMain(object):
    def setupUi(self, CascadesMain):
        CascadesMain.setObjectName("CascadesMain")
        CascadesMain.resize(1167, 622)
        self.verticalLayout = QtWidgets.QVBoxLayout(CascadesMain)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_2 = QtWidgets.QLabel(CascadesMain)
        self.label_2.setStyleSheet("font-size: 24px;\n"
"font-weight: bold;")
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.label = QtWidgets.QLabel(CascadesMain)
        self.label.setStyleSheet("font-weight: bold;")
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.loadExistingButton = QtWidgets.QPushButton(CascadesMain)
        self.loadExistingButton.setObjectName("loadExistingButton")
        self.horizontalLayout_3.addWidget(self.loadExistingButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.label_3 = QtWidgets.QLabel(CascadesMain)
        self.label_3.setStyleSheet("font-weight: bold;")
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_5 = QtWidgets.QLabel(CascadesMain)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_2.addWidget(self.label_5)
        self.selectFastaButton = QtWidgets.QPushButton(CascadesMain)
        self.selectFastaButton.setObjectName("selectFastaButton")
        self.horizontalLayout_2.addWidget(self.selectFastaButton)
        self.selectedFileLabel = QtWidgets.QLabel(CascadesMain)
        self.selectedFileLabel.setMinimumSize(QtCore.QSize(200, 0))
        self.selectedFileLabel.setText("")
        self.selectedFileLabel.setObjectName("selectedFileLabel")
        self.horizontalLayout_2.addWidget(self.selectedFileLabel)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.createCascadeTable = QtWidgets.QTableView(CascadesMain)
        self.createCascadeTable.setObjectName("createCascadeTable")
        self.verticalLayout.addWidget(self.createCascadeTable)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.label_8 = QtWidgets.QLabel(CascadesMain)
        self.label_8.setObjectName("label_8")
        self.horizontalLayout.addWidget(self.label_8)
        self.domainCombo = QtWidgets.QComboBox(CascadesMain)
        self.domainCombo.setMinimumSize(QtCore.QSize(200, 0))
        self.domainCombo.setObjectName("domainCombo")
        self.horizontalLayout.addWidget(self.domainCombo)
        self.label_7 = QtWidgets.QLabel(CascadesMain)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy)
        self.label_7.setObjectName("label_7")
        self.horizontalLayout.addWidget(self.label_7)
        self.stepsSpinBox = QtWidgets.QSpinBox(CascadesMain)
        self.stepsSpinBox.setMinimum(2)
        self.stepsSpinBox.setMaximum(10)
        self.stepsSpinBox.setProperty("value", 4)
        self.stepsSpinBox.setObjectName("stepsSpinBox")
        self.horizontalLayout.addWidget(self.stepsSpinBox)
        self.label_6 = QtWidgets.QLabel(CascadesMain)
        self.label_6.setMinimumSize(QtCore.QSize(0, 0))
        self.label_6.setObjectName("label_6")
        self.horizontalLayout.addWidget(self.label_6)
        self.identityTargetSpinBox = QtWidgets.QSpinBox(CascadesMain)
        self.identityTargetSpinBox.setMinimum(10)
        self.identityTargetSpinBox.setProperty("value", 50)
        self.identityTargetSpinBox.setObjectName("identityTargetSpinBox")
        self.horizontalLayout.addWidget(self.identityTargetSpinBox)
        self.label_4 = QtWidgets.QLabel(CascadesMain)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout.addWidget(self.label_4)
        self.emailLineEdit = QtWidgets.QLineEdit(CascadesMain)
        self.emailLineEdit.setMinimumSize(QtCore.QSize(200, 0))
        self.emailLineEdit.setObjectName("emailLineEdit")
        self.horizontalLayout.addWidget(self.emailLineEdit)
        self.createButton = QtWidgets.QPushButton(CascadesMain)
        self.createButton.setObjectName("createButton")
        self.horizontalLayout.addWidget(self.createButton)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(CascadesMain)
        QtCore.QMetaObject.connectSlotsByName(CascadesMain)

    def retranslateUi(self, CascadesMain):
        _translate = QtCore.QCoreApplication.translate
        CascadesMain.setWindowTitle(_translate("CascadesMain", "Form"))
        self.label_2.setText(_translate("CascadesMain", "Cascade BLAST"))
        self.label.setText(_translate("CascadesMain", "Option 1: Load existing results"))
        self.loadExistingButton.setText(_translate("CascadesMain", "Select File"))
        self.label_3.setText(_translate("CascadesMain", "Option 2: Create new"))
        self.label_5.setText(_translate("CascadesMain", "Select Fasta File"))
        self.selectFastaButton.setText(_translate("CascadesMain", "Select File"))
        self.label_8.setText(_translate("CascadesMain", "Domain"))
        self.label_7.setText(_translate("CascadesMain", "Number of steps in Cascade"))
        self.label_6.setText(_translate("CascadesMain", "minimum identity target"))
        self.label_4.setText(_translate("CascadesMain", "email"))
        self.createButton.setText(_translate("CascadesMain", "Create"))
