# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/msa/visual/ScoreByDivergence.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_ScoreByDivergence(object):
    def setupUi(self, ScoreByDivergence):
        ScoreByDivergence.setObjectName("ScoreByDivergence")
        ScoreByDivergence.resize(981, 620)
        self.verticalLayout = QtWidgets.QVBoxLayout(ScoreByDivergence)
        self.verticalLayout.setObjectName("verticalLayout")
        self.descriptionLabel = QtWidgets.QLabel(ScoreByDivergence)
        self.descriptionLabel.setText("")
        self.descriptionLabel.setWordWrap(True)
        self.descriptionLabel.setObjectName("descriptionLabel")
        self.verticalLayout.addWidget(self.descriptionLabel)
        spacerItem = QtWidgets.QSpacerItem(20, 159, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)

        self.retranslateUi(ScoreByDivergence)
        QtCore.QMetaObject.connectSlotsByName(ScoreByDivergence)

    def retranslateUi(self, ScoreByDivergence):
        _translate = QtCore.QCoreApplication.translate
        ScoreByDivergence.setWindowTitle(_translate("ScoreByDivergence", "Form"))
