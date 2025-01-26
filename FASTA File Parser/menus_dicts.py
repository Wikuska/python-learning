from PySide6.QtGui import QAction

import_menu_actions_dict = {
    "import_file_action": "&Import file from your device"
}

file_menu_actions_dict = {
    "records_number_action": "&Get number of records",
    "records_length_action": "&Get length of each record",
    "shortest_longest_record_action": "&Get the longest and the sortest record"
}

dna_menu_actions_dict = {
    "sequence_length_action": "&Get sequence length",
    "nucleotides_count_action": "&Get number of each nucleotide",
    "gc_content_action": "&Get GC %",
    "complementary_action": "&Get complementary sequence",
    "reverse_complementary_action": "&Get reverse complementary sequence",
    "transcribe_action": "&Transcribe to mRNA"
}

def dict_to_object_dict(app, function, menu):
    objects = {}

    if menu == "import":
        dictionary = import_menu_actions_dict
    elif menu == "file":
        dictionary = file_menu_actions_dict
    elif menu == "dna":
        dictionary = dna_menu_actions_dict

    for action_name, action_descr in dictionary.items():
        action = QAction(action_descr, app)
        action.setData(action_name.replace("_action", ""))
        action.triggered.connect(function)
        objects[action_name] = action

    return objects
