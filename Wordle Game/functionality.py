
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
