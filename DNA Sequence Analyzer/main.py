from flask import Flask, render_template, request

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    if request.method == "POST":
        # Get the action (which button was clicked)
        action = request.form.get("action")
        # Get DNA sequence from the form 
        dna_sequence = request.form["dna_sequence"].strip().upper()
        
        if action == "import":
            result = {"error": "Feature coming soon"}
        elif action == "analyze":
            if not all(base in "ATGC" for base in dna_sequence):
                result = {"error": "Invalid sequence! Use only A, T, G, and C."}
            elif dna_sequence == "":
                result = {"error": "No sequence in form."}
            else:
                base_counts = {
                    "A": dna_sequence.count("A"),
                    "T": dna_sequence.count("T"),
                    "G": dna_sequence.count("G"),
                    "C": dna_sequence.count("C")
                }
                gc_content = ((base_counts["G"] + base_counts["C"]) / len(dna_sequence)) * 100
                result = {
                    "base_counts": base_counts,
                    "gc_content": f"{gc_content:.2f}%"
                }
    
    return render_template("index.html", result=result)

if __name__ == "__main__":
    app.run(debug=True)
