# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'interface5.5.ui'
#
# Created by: PyQt5 UI code generator 5.11.3
#
# WARNING! All changes made in this file will be lost!
"""
This is the interface generated by Qt Designer and converted to a python script
"""

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(628, 626)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.gridLayout.addLayout(self.verticalLayout_4, 4, 0, 1, 1)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.InterrogationAreas = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setItalic(True)
        self.InterrogationAreas.setFont(font)
        self.InterrogationAreas.setObjectName("InterrogationAreas")
        self.verticalLayout_3.addWidget(self.InterrogationAreas)
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_3.addWidget(self.label_3)
        self.WindowSize = QtWidgets.QLineEdit(self.centralwidget)
        self.WindowSize.setObjectName("WindowSize")
        self.verticalLayout_3.addWidget(self.WindowSize)
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_3.addWidget(self.label_4)
        self.overlap = QtWidgets.QLineEdit(self.centralwidget)
        self.overlap.setObjectName("overlap")
        self.verticalLayout_3.addWidget(self.overlap)
        self.gridLayout.addLayout(self.verticalLayout_3, 3, 0, 1, 1)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.Method = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setItalic(True)
        self.Method.setFont(font)
        self.Method.setObjectName("Method")
        self.verticalLayout_2.addWidget(self.Method)
        self.write_method = QtWidgets.QComboBox(self.centralwidget)
        self.write_method.setObjectName("write_method")
        self.write_method.addItem("")
        self.write_method.addItem("")
        self.verticalLayout_2.addWidget(self.write_method)
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setObjectName("label_8")
        self.verticalLayout_2.addWidget(self.label_8)
        self.number_iterations2 = QtWidgets.QSpinBox(self.centralwidget)
        self.number_iterations2.setObjectName("number_iterations2")
        self.verticalLayout_2.addWidget(self.number_iterations2)
        self.gridLayout.addLayout(self.verticalLayout_2, 2, 1, 2, 3)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.Images = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(12)
        font.setItalic(True)
        self.Images.setFont(font)
        self.Images.setObjectName("Images")
        self.verticalLayout.addWidget(self.Images)
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.image1 = QtWidgets.QLineEdit(self.centralwidget)
        self.image1.setObjectName("image1")
        self.verticalLayout.addWidget(self.image1)
        self.image1_button = QtWidgets.QToolButton(self.centralwidget)
        self.image1_button.setObjectName("image1_button")
        self.verticalLayout.addWidget(self.image1_button)
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.Number_of_images = QtWidgets.QSpinBox(self.centralwidget)
        self.Number_of_images.setMaximum(100000000)
        self.Number_of_images.setObjectName("Number_of_images")
        self.verticalLayout.addWidget(self.Number_of_images)
        self.gridLayout.addLayout(self.verticalLayout, 2, 0, 1, 1)
        self.Ok_button2 = QtWidgets.QPushButton(self.centralwidget)
        self.Ok_button2.setObjectName("Ok_button2")
        self.gridLayout.addWidget(self.Ok_button2, 4, 3, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        font = QtGui.QFont()
        font.setPointSize(14)
        font.setItalic(True)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.gridLayout.addWidget(self.label_10, 0, 0, 1, 1)
        self.quit_button = QtWidgets.QPushButton(self.centralwidget)
        self.quit_button.setObjectName("quit_button")
        self.gridLayout.addWidget(self.quit_button, 4, 2, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 628, 29))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionClose = QtWidgets.QAction(MainWindow)
        self.actionClose.setObjectName("actionClose")
        self.actionLoad_Images = QtWidgets.QAction(MainWindow)
        self.actionLoad_Images.setObjectName("actionLoad_Images")
        self.actionUML_Diagram = QtWidgets.QAction(MainWindow)
        self.actionUML_Diagram.setObjectName("actionUML_Diagram")
        self.menuFile.addAction(self.actionClose)
        self.menuFile.addAction(self.actionLoad_Images)
        self.menuHelp.addAction(self.actionUML_Diagram)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.InterrogationAreas.setText(_translate("MainWindow", "Intorrogation Areas"))
        self.label_3.setText(_translate("MainWindow", "Window Size"))
        self.WindowSize.setPlaceholderText(_translate("MainWindow", "Type a value between 16 and image size"))
        self.label_4.setText(_translate("MainWindow", "Overlap"))
        self.overlap.setPlaceholderText(_translate("MainWindow", "Type a value between 16 and image size "))
        self.Method.setText(_translate("MainWindow", "Method"))
        self.write_method.setItemText(0, _translate("MainWindow", "Multigrid"))
        self.write_method.setItemText(1, _translate("MainWindow", "Normal"))
        self.label_8.setText(_translate("MainWindow", "Number of iterations"))
        self.Images.setText(_translate("MainWindow", "Images"))
        self.label.setText(_translate("MainWindow", "Image 1"))
        self.image1_button.setText(_translate("MainWindow", "..."))
        self.label_2.setText(_translate("MainWindow", "Number of images to load"))
        self.Ok_button2.setText(_translate("MainWindow", "Ok"))
        self.label_10.setText(_translate("MainWindow", "PIV Method Implementation Code"))
        self.quit_button.setText(_translate("MainWindow", "Quit"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionLoad_Images.setText(_translate("MainWindow", "Load Images"))
        self.actionUML_Diagram.setText(_translate("MainWindow", "UML Diagram"))

