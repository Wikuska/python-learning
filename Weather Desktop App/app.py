from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMainWindow, QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QLabel, QComboBox
from funtionality import get_weather, get_icon

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Weather App")
        self.setFixedSize(380, 420)
        self.setStyleSheet("background-color: #a2d2ff")

        layout = QVBoxLayout()
        layout.setContentsMargins(10,20,10,0)

        top_layout = QVBoxLayout()
        output_layout = QVBoxLayout()
        output_layout.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
        output_layout.setSpacing(10)

        main_label = QLabel("Enter city name", self, alignment=Qt.AlignmentFlag.AlignHCenter)
        main_label.setStyleSheet("""
            font-size: 20px;
            font-weight: bold;
            font-family: 'Arial', sans-serif;
        """)

        self.city_input = QLineEdit(self, alignment=Qt.AlignmentFlag.AlignHCenter)
        self.city_input.setPlaceholderText("Enter city name")
        self.city_input.setStyleSheet("background-color: white")
        self.city_input.setFixedHeight(30)

        self.get_weather_button = QPushButton("Get Weather", self)
        self.get_weather_button.setFixedSize(150, 30)
        self.get_weather_button.setStyleSheet("background-color: white")
        self.get_weather_button.clicked.connect(self.get_weather)

        self.city_label = QLabel(self)
        self.city_label.setStyleSheet("""
            font-size: 20px;
            font-weight: bold;
            font-family: 'Arial', sans-serif;
        """)
        self.icon_label = QLabel(self)
        self.temp_label = QLabel(self)
        self.temp_label.setStyleSheet("font-size: 13px;")
        self.weather_label = QLabel(self)
        self.weather_label.setStyleSheet("font-size: 13px;")

        top_layout.addWidget(main_label)
        top_layout.addSpacing(10)
        top_layout.addWidget(self.city_input)
        top_layout.addWidget(self.get_weather_button, alignment=Qt.AlignmentFlag.AlignHCenter)
        output_layout.addWidget(self.city_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        output_layout.addWidget(self.icon_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        output_layout.addWidget(self.temp_label, alignment=Qt.AlignmentFlag.AlignHCenter)
        output_layout.addWidget(self.weather_label, alignment=Qt.AlignmentFlag.AlignHCenter)

        layout.addLayout(top_layout, stretch=0)
        layout.setAlignment(top_layout, Qt.AlignmentFlag.AlignTop)
        layout.addLayout(output_layout, stretch=1)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def get_weather(self):
        if self.city_input.text() == "": 
            return
        weather_id, weather, temp_in_c = get_weather(self.city_input.text())
        self.city_label.setText(f"Weather in {self.city_input.text().title()}")
        self.icon_label.setPixmap(get_icon(weather_id))
        self.temp_label.setText(f"Temperature: {temp_in_c}°C")
        self.weather_label.setText(f"Weather: {weather.capitalize()}")


app = QApplication([])
w = MainWindow()
w.show()
app.exec()
