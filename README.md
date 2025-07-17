# AI-Enhanced Peptide Generator

A web-based app to upload protein structures (.pdb), extract sequence and spatial info, and generate AI-suggested peptide candidates with explanations.

## Features
- Upload PDB files
- Extract sequence and residue coordinates
- Generate and explain peptide candidates using LLMs
- Modular, beginner-friendly codebase

## Setup
1. Clone this repo
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the app:
   ```bash
   streamlit run app.py
   ```

## Usage
- Upload a `.pdb` file
- Select chain and parameters
- Enter your LLM API key (OpenAI, Claude, etc.)
- View extracted info and peptide suggestions

## Tech Stack
- Python 3.10+
- Streamlit UI
- Biopython for structure parsing
- OpenAI/LLM API support

---

*Each module is testable and beginner-friendly. See code comments for details.* 