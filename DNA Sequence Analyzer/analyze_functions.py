
def validate_sequence(sequence):
    if not all(base in "ATGC" for base in sequence) or sequence == "":
        return False
    else:
        return True
    
def base_analyze(sequence):
    base_counts = {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }
    gc_content = ((base_counts["G"] + base_counts["C"]) / len(sequence)) * 100
    result = f"A: {base_counts["A"]}\nT: {base_counts["T"]}\nG: {base_counts["G"]}\nC: {base_counts["C"]}\nGC Content: {gc_content:.2f}%"
    return result

def rna_transcription(sequence):
    transcription = sequence.replace("T", "U")
    return transcription
