from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout
from Funcs import create_board
from PlayerNameDialog import PlayerNameDialog


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Battleship Game")
        
        dialog = PlayerNameDialog()
        if dialog.exec():
            p1_name = dialog.p1_input.text() or dialog.p1_input.placeholderText()
            p2_name = dialog.p2_input.text() or dialog.p2_input.placeholderText()

        main_layout = QVBoxLayout()
        boards_layout = QHBoxLayout()
        
        p1_board, p1_buttons = create_board(p1_name)
        p2_board, p2_buttons = create_board(p2_name)
        
        boards_layout.addLayout(p1_board)
        boards_layout.addSpacing(50)
        boards_layout.addLayout(p2_board)
        
        main_layout.addLayout(boards_layout)
        
        print(p1_name)
        print(p2_name)
        
        self.setLayout(main_layout)
        
app = QApplication([])
w = MainWindow()
w.show()
app.exec()       