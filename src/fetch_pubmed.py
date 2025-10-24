"""
PubMed Fetcher Script
---------------------

A universal Python tool to fetch and parse PubMed articles using Biopython and NCBI Entrez API.
Search for any topic, extract titles, authors, abstracts, and export results to CSV.

Author: Aarti Nagpal
Email: aarti96028@gmail.com
Date: October 2025
"""

from Bio import Entrez
import pandas as pd
from pathlib import Path
from datetime import datetime
import time

def fetch_pubmed(term: str, email: str, retmax: int = 10):
    """Fetch PubMed abstracts for a given search term."""
    
    Entrez.email = email
    Entrez.tool = "PubMedFetcher"
    print(f"\n🔍 Searching PubMed for: '{term}'")

    # Search PubMed
    search_results = Entrez.esearch(db="pubmed", term=term, retmax=retmax)
    search_results = Entrez.read(search_results)
    total_count = int(search_results["Count"])
    print(f"📊 Found {total_count} total results. Fetching top {len(search_results['IdList'])}...")

    # Fetch the records in XML format
    records = Entrez.efetch(
        db="pubmed",
        id=search_results['IdList'],
        rettype="abstract",
        retmode="xml"
    )
    parsed_records = Entrez.read(records)
    records.close()

    articles_list = []

    # Iterate through fetched articles
    for article in parsed_records["PubmedArticle"]:
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

            # Create dictionary for each article
            article_dict = {
                'PMID': medline_citation.get('PMID', 'N/A'),
                'Title': article_info.get('ArticleTitle', 'No Title'),
                'Authors': '; '.join(authors) if authors else 'No Authors',
                'Abstract': abstract_text
            }
            articles_list.append(article_dict)
        
        except Exception as e:
            print(f"⚠️ Error processing article: {e}")

        time.sleep(0.3)  # Pause to comply with NCBI rate limits

    # Create a DataFrame
    df = pd.DataFrame(articles_list)

    # Save the results into a CSV file with timestamp
    Path("data").mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    output_path = f"data/pubmed_results_{timestamp}.csv"
    df.to_csv(output_path, index=False, encoding='utf-8')

    print(f"\n✅ {len(df)} articles fetched successfully!")
    print(f"💾 Results saved to: {output_path}")

    return df


if __name__ == "__main__":
    # Example usage
    df = fetch_pubmed(term="Machine learning in biology", email="aarti96028@gmail.com", retmax=10)
    print("\n", df.head())
