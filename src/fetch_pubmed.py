
---

```python
"""
fetch_pubmed.py

A universal Python script to automatically fetch and parse PubMed articles
using the NCBI Entrez API via Biopython.

This tool allows users to search PubMed using any query term and extract
metadata such as Title, Authors, Abstract, and PMID into a structured
Pandas DataFrame or CSV file.

Author: Aarti Nagpal
Email: aarti96028@gmail.com
Date: October 2025
"""

from Bio import Entrez
import pandas as pd
import argparse
import sys
from tqdm import tqdm

def fetch_pubmed_articles(term, email, retmax=20, save_csv=True):
    """Fetch PubMed articles based on a search term."""

    Entrez.email = email
    print(f"🔍 Searching PubMed for: '{term}'")

    try:
        # Search PubMed
        search_handle = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        ids = search_results.get("IdList", [])
        print(f"📄 Found {len(ids)} articles.")

        if not ids:
            print("No articles found for the given query.")
            return pd.DataFrame()

        # Fetch article details
        fetch_handle = Entrez.efetch(db="pubmed", id=ids, rettype="abstract", retmode="xml")
        records = Entrez.read(fetch_handle)
        fetch_handle.close()

        articles_list = []

        for article in tqdm(records["PubmedArticle"], desc="Processing articles"):
            try:
                medline_citation = article.get('MedlineCitation', {})
                article_info = medline_citation.get('Article', {})

                # Extract Authors
                authors = []
                if 'AuthorList' in article_info:
                    for author in article_info['AuthorList']:
                        if isinstance(author, dict) and 'LastName' in author and 'ForeName' in author:
                            authors.append(f"{author['LastName']} {author['ForeName']}")

                # Extract Abstract
                abstract_text = " ".join(article_info.get('Abstract', {}).get('AbstractText', ["No Abstract"]))

                article_dict = {
                    'PMID': medline_citation.get('PMID', 'N/A'),
                    'Title': article_info.get('ArticleTitle', 'No Title'),
                    'Authors': '; '.join(authors) if authors else 'No Authors',
                    'Abstract': abstract_text
                }
                articles_list.append(article_dict)
            
            except Exception as e:
                print(f"Error processing article: {e}")

        df = pd.DataFrame(articles_list)
        print("\n Articles successfully fetched!")

        if save_csv:
            output_path = "data/pubmed_results.csv"
            df.to_csv(output_path, index=False)
            print(f" Saved results to: {output_path}")

        return df

    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch PubMed abstracts for any biomedical topic.")
    parser.add_argument("--term", type=str, required=True,
                        help="Search term for PubMed query (e.g., 'cancer immunotherapy').")
    parser.add_argument("--email", type=str, required=True, help="Your email address (required by NCBI).")
    parser.add_argument("--retmax", type=int, default=20, help="Number of articles to retrieve.")
    parser.add_argument("--no-csv", action="store_true", help="Do not save output to CSV file.")

    args = parser.parse_args()

    fetch_pubmed_articles(args.term, args.email, args.retmax, not args.no_csv)
