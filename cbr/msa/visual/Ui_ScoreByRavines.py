# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/msa/visual/ScoreByRavines.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_ScoreByRavines(object):
    def setupUi(self, ScoreByRavines):
        ScoreByRavines.setObjectName("ScoreByRavines")
        ScoreByRavines.resize(400, 300)
        self.gridLayout = QtWidgets.QGridLayout(ScoreByRavines)
        self.gridLayout.setObjectName("gridLayout")
        self.continuityTresholdLabel = QtWidgets.QLabel(ScoreByRavines)
        self.continuityTresholdLabel.setMinimumSize(QtCore.QSize(50, 0))
        self.continuityTresholdLabel.setMaximumSize(QtCore.QSize(50, 16777215))
        self.continuityTresholdLabel.setText("")
        self.continuityTresholdLabel.setObjectName("continuityTresholdLabel")
        self.gridLayout.addWidget(self.continuityTresholdLabel, 1, 2, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 237, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 2, 1, 1, 1)
        self.continuityTresholdSlider = QtWidgets.QSlider(ScoreByRavines)
        self.continuityTresholdSlider.setProperty("value", 80)
        self.continuityTresholdSlider.setOrientation(QtCore.Qt.Horizontal)
        self.continuityTresholdSlider.setObjectName("continuityTresholdSlider")
        self.gridLayout.addWidget(self.continuityTresholdSlider, 1, 0, 1, 2)
        self.descriptionLabel = QtWidgets.QLabel(ScoreByRavines)
        self.descriptionLabel.setText("")
        self.descriptionLabel.setWordWrap(True)
        self.descriptionLabel.setObjectName("descriptionLabel")
        self.gridLayout.addWidget(self.descriptionLabel, 0, 0, 1, 3)

        self.retranslateUi(ScoreByRavines)
        QtCore.QMetaObject.connectSlotsByName(ScoreByRavines)

    def retranslateUi(self, ScoreByRavines):
        _translate = QtCore.QCoreApplication.translate
        ScoreByRavines.setWindowTitle(_translate("ScoreByRavines", "Form"))
