
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
    return True
