from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMainWindow, QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QLabel, QComboBox

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Weather App")
        self.setFixedSize(380, 480)

        layout = QVBoxLayout()
        layout.setContentsMargins(20,20,20,20)

        top_layout = QVBoxLayout()
        output_layout = QVBoxLayout()

        self.city_input = QLineEdit(self)
        self.city_input.setPlaceholderText("Enter city name")
        self.city_input.setFixedHeight(30)

        self.get_weather_button = QPushButton("Get Weather", self)
        self.get_weather_button.setFixedSize(150, 30)

        self.icon_label = QLabel(self)
        self.temp_label = QLabel(self)
        self.weather_label = QLabel(self)

        top_layout.addWidget(self.city_input)
        top_layout.addWidget(self.get_weather_button, alignment=Qt.AlignmentFlag.AlignHCenter)

        output_layout.addWidget(self.icon_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        output_layout.addWidget(self.temp_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        output_layout.addWidget(self.weather_label, alignment=Qt.AlignmentFlag.AlignHCenter)

        layout.addLayout(top_layout)
        layout.setAlignment(top_layout, Qt.AlignmentFlag.AlignTop)
        layout.addSpacing(50)
        layout.addLayout(output_layout)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

app = QApplication([])
w = MainWindow()
w.show()
app.exec()
