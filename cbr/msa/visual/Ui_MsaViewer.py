# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MsaViewer.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MsaViewer(object):
    def setupUi(self, MsaViewer):
        MsaViewer.setObjectName("MsaViewer")
        MsaViewer.resize(400, 300)
        self.verticalLayout = QtWidgets.QVBoxLayout(MsaViewer)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_2 = QtWidgets.QLabel(MsaViewer)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.structuresCombo = QtWidgets.QComboBox(MsaViewer)
        self.structuresCombo.setObjectName("structuresCombo")
        self.horizontalLayout.addWidget(self.structuresCombo)
        self.structuresRefreshButton = QtWidgets.QPushButton(MsaViewer)
        self.structuresRefreshButton.setObjectName("structuresRefreshButton")
        self.horizontalLayout.addWidget(self.structuresRefreshButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(MsaViewer)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        self.selectMsaButton = QtWidgets.QPushButton(MsaViewer)
        self.selectMsaButton.setObjectName("selectMsaButton")
        self.horizontalLayout_2.addWidget(self.selectMsaButton)
        self.createMsaButton = QtWidgets.QPushButton(MsaViewer)
        self.createMsaButton.setObjectName("createMsaButton")
        self.horizontalLayout_2.addWidget(self.createMsaButton)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_3 = QtWidgets.QLabel(MsaViewer)
        self.label_3.setMinimumSize(QtCore.QSize(136, 0))
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_3.addWidget(self.label_3)
        self.sequenceCombo = QtWidgets.QComboBox(MsaViewer)
        self.sequenceCombo.setObjectName("sequenceCombo")
        self.horizontalLayout_3.addWidget(self.sequenceCombo)
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem2)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.colorConservedButton = QtWidgets.QPushButton(MsaViewer)
        self.colorConservedButton.setMinimumSize(QtCore.QSize(151, 0))
        self.colorConservedButton.setObjectName("colorConservedButton")
        self.horizontalLayout_4.addWidget(self.colorConservedButton)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem3)
        self.verticalLayout.addLayout(self.horizontalLayout_4)
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem4)

        self.retranslateUi(MsaViewer)
        QtCore.QMetaObject.connectSlotsByName(MsaViewer)

    def retranslateUi(self, MsaViewer):
        _translate = QtCore.QCoreApplication.translate
        MsaViewer.setWindowTitle(_translate("MsaViewer", "Form"))
        self.label_2.setText(_translate("MsaViewer", "Structure"))
        self.structuresRefreshButton.setText(_translate("MsaViewer", "Refresh"))
        self.label.setText(_translate("MsaViewer", "MSA File"))
        self.selectMsaButton.setText(_translate("MsaViewer", "Browse"))
        self.createMsaButton.setText(_translate("MsaViewer", "Create"))
        self.label_3.setText(_translate("MsaViewer", "MSA Structure Sequence"))
        self.colorConservedButton.setText(_translate("MsaViewer", "Color Conserved Regions"))
