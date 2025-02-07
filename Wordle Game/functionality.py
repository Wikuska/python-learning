from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QHBoxLayout, QVBoxLayout
import requests
import random

def get_word():
    response = requests.get("https://api.datamuse.com/words?ml=common&max=100")
    words_data = [word["word"] for word in response.json()]
    filtered_words = [word for word in words_data if len(word) == 5 and len(word.split()) == 1]
    random_word = random.choice(filtered_words)
    return random_word

def create_labels(main_layout):
    label_obj_list = []
    rows_layout = QVBoxLayout()
    for _ in range(6):
        row_layout = QHBoxLayout()
        for _ in range(5):
            label = QLabel("", alignment = Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter)
            label.setStyleSheet("font-size: 24px; font-weight: bold; background-color: #3a3a3c")
            label_obj_list.append(label)
            row_layout.addWidget(label)
        rows_layout.addLayout(row_layout)

    if main_layout.count() > 0:
        item = main_layout.takeAt(0)
        if item:
            layout = item.layout()
            if layout:
                while layout.count():
                    child = layout.takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
                layout.deleteLater()

    main_layout.insertLayout(0, rows_layout)

    return label_obj_list

def validate_answer(guess, reason_label):
    if guess == "":
        reason_label.setText("No guess given")
        return False
    if len(guess.split()) > 1:
        reason_label.setText("Use one word only")
        return False
    if len(guess) != 5:
        reason_label.setText("Use 5 letter word")
        return False
    reason_label.setText("")
    return True

def check_correct_letters(guess, answer, labels_list):

    letters_score = [] # wrong letter - 0, good letter wrong placement - 1, good letter good placement - 2

    for gletter, aletter in zip(guess, answer):  
        if gletter in answer:
            if gletter == aletter:
                letters_score.append(2)
            else:
                letters_score.append(1)
        else:
            letters_score.append(0)

    for label, letter, score in zip(labels_list, guess, letters_score):
        label.setText(letter)
        if score == 1:
            label.setStyleSheet("font-size: 24px; font-weight: bold; background-color: yellow")
        elif score == 2:
            label.setStyleSheet("font-size: 24px; font-weight: bold; background-color: green")
    del labels_list[:5]

    if guess == answer:
        return True
    else:
        return False
