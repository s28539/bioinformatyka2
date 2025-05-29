#!/usr/bin/env python3
"""
Enhanced NCBI GenBank Data Retriever
- Filtrowanie sekwencji wg długości
- Generowanie CSV
- Tworzenie wykresu długości sekwencji
"""

from Bio import Entrez, SeqIO
import time
import csv
import matplotlib.pyplot as plt
from Bio.Entrez import api_key


class NCBIRetriever:
    def __init__(self, email, api_key):
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'EnhancedRetriever'

    def search_taxid(self, taxid, min_len, max_len):
        print(f"🔍 Searching for records with taxID: {taxid} and length {min_len}-{max_len} bp")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"🧬 Organism: {organism_name} (TaxID: {taxid})")

            # ✅ Search term with sequence length filter
            search_term = f"txid{taxid}[Organism] AND {min_len}:{max_len}[Sequence Length]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])

            if count == 0:
                print("⚠️ No records found.")
                return None

            print(f"✅ Found {count} records")

            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count
            return count

        except Exception as e:
            print(f"❌ Error searching TaxID {taxid}: {e}")
            return None

    def fetch_records(self, max_records=100):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("⚠️ No search results. Run search_taxid() first.")
            return []

        records = []
        batch_size = 500
        for start in range(0, min(self.count, max_records), batch_size):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=start,
                    retmax=min(batch_size, max_records - start),
                    webenv=self.webenv,
                    query_key=self.query_key
                )
                batch_records = list(SeqIO.parse(handle, "gb"))
                records.extend(batch_records)
                print(f"📥 Fetched {len(batch_records)} records (batch starting at {start})")
                time.sleep(0.5)
            except Exception as e:
                print(f"❌ Error fetching batch at {start}: {e}")
                break
        return records

    def filter_records_by_length(self, records, min_len, max_len):
        filtered = [r for r in records if min_len <= len(r.seq) <= max_len]
        print(f"🔎 Filtered {len(filtered)} records between {min_len} and {max_len} bp")
        return filtered

    def save_to_csv(self, records, filename):
        with open(filename, mode='w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Accession", "Length", "Description"])
            for record in records:
                writer.writerow([record.id, len(record.seq), record.description])
        print(f"💾 Saved CSV to {filename}")

    def plot_lengths(self, records, filename):
        if not records:
            print("⚠️ No records to plot.")
            return

        sorted_records = sorted(records, key=lambda r: len(r.seq), reverse=True)
        accessions = [r.id for r in sorted_records]
        lengths = [len(r.seq) for r in sorted_records]

        plt.figure(figsize=(max(12, len(accessions) * 0.5), 6))  # dynamic width

        # Wykres liniowy + punkty
        x_positions = list(range(len(accessions)))
        plt.plot(x_positions, lengths, linestyle='-', color='skyblue', label='Length')
        plt.scatter(x_positions, lengths, color='blue', s=40, zorder=3, label='Points')

        # Etykiety X: wszystkie accession IDs
        plt.xticks(x_positions, accessions, rotation=90, fontsize=6)

        plt.xlabel("GenBank Accession Number")
        plt.ylabel("Sequence Length (bp)")
        plt.title("Sequence Lengths by Accession (sorted)")
        plt.grid(True, linestyle='--', alpha=0.4)

        # Opcjonalnie: długości nad punktami jeśli <= 25 rekordów
        if len(records) <= 25:
            for i, length in enumerate(lengths):
                plt.text(i, length + 20, str(length), ha='center', va='bottom', fontsize=7, rotation=45)

        plt.tight_layout()
        plt.savefig(filename)
        plt.close()
        print(f"📊 Saved plot to {filename}")


def main():
    print("🧪 NCBI GenBank Sequence Retriever\n")
    # email = input("Enter your email address for NCBI: ")
    # api_key = input("Enter your NCBI API key: ")
    # taxid = input("Enter taxonomic ID (taxid): ")
    email = "s28539@pjwstk.edu.pl"
    api_key = "ba7a640fb0abc6bb839be42d04971e380608"
    taxid = 9606

    try:
        min_len = int(input("Minimum sequence length (bp): "))
        max_len = int(input("Maximum sequence length (bp): "))
        max_records = int(input("Max number of records to download: "))
    except ValueError:
        print("❌ Invalid input. Please enter numeric values.")
        return

    retriever = NCBIRetriever(email, api_key)

    if not retriever.search_taxid(taxid, min_len, max_len):
        return

    print("\n📡 Fetching records from NCBI...")
    all_records = retriever.fetch_records(max_records=max_records)

    print("\n🧮 Filtering records by sequence length...")
    filtered = retriever.filter_records_by_length(all_records, min_len, max_len)

    if not filtered:
        print("⚠️ No records matched the length criteria.")
        return

    base_filename = f"taxid_{taxid}_filtered"

    retriever.save_to_csv(filtered, f"{base_filename}.csv")
    retriever.plot_lengths(filtered, f"{base_filename}_plot.png")

    print("\n✅ Done!")


if __name__ == "__main__":
    main()
