# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Main.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Main(object):
    def setupUi(self, Main):
        Main.setObjectName("Main")
        Main.resize(767, 656)
        Main.setStyleSheet("")
        self.verticalLayout = QtWidgets.QVBoxLayout(Main)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.widget = QtWidgets.QWidget(Main)
        self.widget.setMinimumSize(QtCore.QSize(175, 130))
        self.widget.setStyleSheet("image: url(:/resources/tumcs.png);")
        self.widget.setObjectName("widget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.widget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.horizontalLayout_2.addWidget(self.widget)
        self.label = QtWidgets.QLabel(Main)
        font = QtGui.QFont()
        font.setPointSize(22)
        font.setBold(True)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.schemaRasppButton = QtWidgets.QPushButton(Main)
        self.schemaRasppButton.setMinimumSize(QtCore.QSize(128, 128))
        font = QtGui.QFont()
        font.setBold(True)
        self.schemaRasppButton.setFont(font)
        self.schemaRasppButton.setStyleSheet("image: url(:/resources/schema.png);")
        self.schemaRasppButton.setIconSize(QtCore.QSize(128, 128))
        self.schemaRasppButton.setFlat(False)
        self.schemaRasppButton.setObjectName("schemaRasppButton")
        self.horizontalLayout_3.addWidget(self.schemaRasppButton)
        self.msaViewerButton = QtWidgets.QPushButton(Main)
        self.msaViewerButton.setMinimumSize(QtCore.QSize(128, 128))
        font = QtGui.QFont()
        font.setBold(True)
        self.msaViewerButton.setFont(font)
        self.msaViewerButton.setStyleSheet("image: url(:/resources/msa.png);")
        self.msaViewerButton.setIconSize(QtCore.QSize(128, 128))
        self.msaViewerButton.setObjectName("msaViewerButton")
        self.horizontalLayout_3.addWidget(self.msaViewerButton)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem2)

        self.retranslateUi(Main)
        QtCore.QMetaObject.connectSlotsByName(Main)

    def retranslateUi(self, Main):
        _translate = QtCore.QCoreApplication.translate
        Main.setWindowTitle(_translate("Main", "Main"))
        self.label.setText(_translate("Main", "CBR Bioinformatics Tools"))
        self.schemaRasppButton.setText(_translate("Main", "SCHEMA-RASPP"))
        self.msaViewerButton.setText(_translate("Main", "MSA Viewer"))
