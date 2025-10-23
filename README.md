
 🧬 **PubMed Fetcher**

A **universal Python tool** to automatically fetch PubMed articles and abstracts using the **NCBI Entrez API** via **Biopython**.  
Use it to retrieve biomedical literature on *any* topic — from cancer research to microbiome studies — and extract metadata for analysis or text mining.

---

🚀 Features
- 🔍 Search PubMed using any keyword or Boolean query  
- 🧠 Extract article **Title**, **Authors**, **Abstract**, and **PMID**  
- 💾 Save results as CSV for downstream analysis  
- 🧩 Integrate easily into bioinformatics or machine learning workflows  

---

💻 Usage

1️⃣Install dependencies

pip install -r requirements.txt

2️⃣ Run the script
python src/fetch_pubmed.py --term "cancer immunotherapy" --email "youremail@example.com" --retmax 10

3️⃣ Output
A CSV file is saved to data/pubmed_results.csv containing:
| PMID     | Title                         | Authors        | Abstract                         |
| -------- | ----------------------------- | -------------- | -------------------------------- |
| 12345678 | Role of cytokines in immunity | Doe J; Smith A | Cytokines play a key role in ... |

🧠 **Example Queries**

"machine learning in bioinformatics"

"CRISPR gene editing"

"microbiome AND metabolic syndrome"

"COVID-19 vaccine response"

🧪 **Applications**

Biomedical literature curation

NLP / BioBERT-based entity extraction

Training datasets for text mining models

Systematic review automation

⚙️ **Requirements**

Python ≥ 3.8

Biopython

Pandas

tqdm

Install them via:
pip install biopython pandas tqdm

🧩 **Project Structure**
pubmed-fetcher/
├── src/
│   └── fetch_pubmed.py
├── data/
│   └── pubmed_results.csv
├── requirements.txt
└── README.md

📄 License
This project is licensed under the MIT License

👩‍🔬 Author
Aarti Nagpal
📧 aarti96028@gmail.com
