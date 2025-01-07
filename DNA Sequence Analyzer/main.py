from flask import Flask, render_template, request
from analyze_functions import validate_sequence, base_analyze, rna_transcription

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    dna_sequence = ""
    if request.method == "POST":
        # Get the action (which button was clicked)
        action = request.form.get("action")
        # Get DNA sequence from the form 
        dna_sequence = request.form["dna_sequence"].strip().upper()
        
        if action == "import":
            result = {"error": "Feature coming soon"}

        elif action == "analyze":
            sequence_validated = validate_sequence(dna_sequence)
            if not sequence_validated:
                result = {"error": "Invalid sequence!"}
            else:
                result = base_analyze(dna_sequence)

        elif action == "to_rna":
            sequence_validated = validate_sequence(dna_sequence)
            if not sequence_validated:
                result = {"error": "Invalid sequence!"}
            else:
                result = rna_transcription(dna_sequence)
    
    return render_template("index.html", result=result, dna_sequence = dna_sequence)

if __name__ == "__main__":
    app.run(debug=True)
