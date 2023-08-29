from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel
from PyQt5.QtGui import QFont, QPixmap
from PyQt5.QtCore import Qt
from typing import Callable

class AppIcon(QWidget):
    def __init__(
            self,
            text : str,
            image_path : str,
            init_app : Callable[[], None]
        ):
        super().__init__()

        self.__init_app = init_app
        self.text = text
        self.image_path = image_path

        self.init_ui()

    def init_ui(self):
        self.setFixedSize(160, 160)
        self.setWindowTitle('Custom Clickable Widget')

        layout = QVBoxLayout()

        pixmap = QPixmap(self.image_path)
        pixmap = pixmap.scaledToWidth(112, Qt.SmoothTransformation)
        image_label = QLabel()
        image_label.setPixmap(pixmap)
        image_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(image_label)

        text_label = QLabel(self.text)
        text_label.setAlignment(Qt.AlignCenter)

        font = QFont()
        font.setPointSize(14)
        font.setBold(True)
        text_label.setFont(font)

        layout.addWidget(text_label)

        self.setLayout(layout)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            print(f"Starting app: {self.text}")
            self.__init_app()
            self.update_clicked(True)

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.update_clicked(False)

    def update_clicked(self, clicked : bool):
        if clicked:
            self.setStyleSheet("background-color: blue;")
        else:
            self.setStyleSheet("")

    def update_hover(self, hovered: bool):

        if hovered:
            self.setStyleSheet("background-color: lightblue;")
        else:
            self.setStyleSheet("")

    def enterEvent(self, event):
        self.update_hover(True)

    def leaveEvent(self, event):
        self.update_hover(False)