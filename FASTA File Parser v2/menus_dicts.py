from PySide6.QtGui import QAction

starting_operations = [
    "Import file - Import FASTA file from device",
    "Work on file - Avaliable only after importing a file. Let you analyse its contetst",
    "Work on sequence - Avaliable only after importing a file. Let you analyse sequence chosen sequence"
]

file_operations = ["Get number of records in file", "Get length of each record", "Get the longest and the sortest record"]
dna_operations = ["Get sequence length", "Get number of each nucleotide", "Get GC %", "Get complementary sequence", "Get reverse complementary sequence", "Transcribe to mRNA"]

main_menu_actions_dict = {
    "import_action": "Import file",
    "file_menu_action": "Work on file",
    "sequence_menu_action": "Work on sequence"
}

def dict_to_object_dict(app, function, menu):
    objects = {}
    if menu == "main":
        dictionary = main_menu_actions_dict

    for action_name, action_descr in dictionary.items():
        action = QAction(action_descr, app)
        action.setData(action_name.replace("_action", ""))
        action.triggered.connect(function)
        objects[action_name] = action

    return objects