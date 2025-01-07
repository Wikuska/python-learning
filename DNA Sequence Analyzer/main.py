from tkinter import *

root = Tk()
root.title("DNA Sequence Analyzer")
root.geometry("400x300")

Label(root, text = "Enter DNA Sequence:").pack(pady = 10)

dna_entry = Entry(root, width = 30)
dna_entry.pack(pady = 5)

def analyze_dna():
    dna_sequence = dna_entry.get().strip().upper()
    if not all(base in "ATGC" for base in dna_sequence):
        print("Invalid sequence! DNA should only contain A, T, G and C.")
        return   

    base_counts = {
        "A": dna_sequence.count("A"),
        "T": dna_sequence.count("T"),
        "G": dna_sequence.count("G"),
        "C": dna_sequence.count("C")
    }

    gc_content = ((base_counts['G'] + base_counts['C']) / len(dna_sequence) * 100)

    result_text = (
        f"Base Counts:\n"
        f"A: {base_counts['A']}\n"
        f"T: {base_counts['T']}\n"
        f"G: {base_counts['G']}\n"
        f"C: {base_counts['C']}\n\n"
        f"GC Content: {gc_content:.2f}%"
    )
    results_label.config(text = result_text)


Button(root, text = "Analyze", command = analyze_dna).pack(pady = 10)

results_label = Label(root, text = "", justify = "left")
results_label.pack(pady = 10)

root.mainloop()
    
