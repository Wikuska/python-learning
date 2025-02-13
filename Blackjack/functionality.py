from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel
import random

def get_hidden_card(card, deck):
    random.shuffle(deck)
    drawn_card = deck.pop()

    card.setText(drawn_card)

    if drawn_card.split()[1] == '♥' or drawn_card.split()[1] == '♦':
        card.setStyleSheet(f"font-size: 20px; color: red; background-color: white")
    else:
        card.setStyleSheet(f"font-size: 20px; color: black; background-color: white")

    return drawn_card.split()[0]

def get_random_card(layout, hidden, deck):
    if hidden:
        card = QLabel(alignment = Qt.AlignmentFlag.AlignCenter)
        card.setFixedSize(100,170)
        card.setStyleSheet(f"font-size: 20px; color: black; background-color: gray")
        card.setFrameShape(QLabel.Box)
        card.setLineWidth(3)

        layout.addWidget(card)
        return card

    random.shuffle(deck)
    drawn_card = deck.pop()

    if drawn_card.split()[1] == '♥' or drawn_card.split()[1] == '♦':
        color = "red"
    else:
        color = "black"

    card = QLabel(drawn_card, alignment = Qt.AlignmentFlag.AlignCenter)
    card.setFixedSize(100,170)
    card.setStyleSheet(f"font-size: 20px; color: {color}; background-color: white")
    card.setFrameShape(QLabel.Box)
    card.setLineWidth(3)

    layout.addWidget(card)

    return drawn_card.split()[0]

def count_score(ranks):
    score = 0

    for rank in ranks:
        try:
            points = int(rank)
        except ValueError:
            if rank == "A":
                points = 11
            else:
                points = 10
        score += points

    if score > 21:
        for rank in ranks:
            if rank == "A":
                score -= 10
                if score <= 21:
                    break

    return score

def clear_layout(layout):
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget:
            widget.deleteLater()
