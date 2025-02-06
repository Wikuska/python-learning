from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMainWindow, QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QLabel, QComboBox
from funtionality import get_weather, get_icon

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Weather App")
        self.setFixedSize(380, 480)
        self.setStyleSheet("background-color: #a2d2ff")

        layout = QVBoxLayout()
        layout.setContentsMargins(20,20,20,20)

        top_layout = QVBoxLayout()
        output_layout = QVBoxLayout()

        main_label = QLabel("Enter city name", self)

        self.city_input = QLineEdit(self, alignment=Qt.AlignmentFlag.AlignHCenter)
        self.city_input.setPlaceholderText("Enter city name")
        self.city_input.setStyleSheet("background-color: white")
        self.city_input.setFixedHeight(30)

        self.get_weather_button = QPushButton("Get Weather", self)
        self.get_weather_button.setFixedSize(150, 30)
        self.get_weather_button.setStyleSheet("background-color: white")
        self.get_weather_button.clicked.connect(self.get_weather)

        self.icon_label = QLabel(self)
        self.temp_label = QLabel(self)
        self.weather_label = QLabel(self)

        top_layout.addWidget(main_label)
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

    def get_weather(self):
        if self.city_input.text() == "": 
            return
        weather_id, weather, temp_in_c = get_weather(self.city_input.text())
        self.icon_label.setPixmap(get_icon(weather_id))
        self.temp_label.setText(f"Temperature: {temp_in_c}°C")
        self.weather_label.setText(f"Weather: {weather.capitalize()}")


app = QApplication([])
w = MainWindow()
w.show()
app.exec()
