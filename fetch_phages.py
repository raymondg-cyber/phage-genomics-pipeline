import os
import time
import ssl
from Bio import Entrez

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context

Entrez.email = "your.email@example.com"

def fetch_unique_klebsiella_phages(target_count=10, output_dir="ncbi_dataset/data"):
    print("📡 Querying NCBI Nucleotide database with universal text-matching parameters...")
    search_term = 'Klebsiella phage AND "complete genome"[Title]'
    
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=60)
        search_results = Entrez.read(handle)
        handle.close()
        
        id_list = search_results["IdList"]
        if not id_list:
            print("❌ No matching phages found.")
            return
            
        print(f"🔍 Found a pool of {len(id_list)} candidate IDs. Filtering for unique entries...")
        os.makedirs(output_dir, exist_ok=True)
        
        downloaded_count = 0
        seen_accessions = set()

        for ncb_id in id_list:
            if downloaded_count >= target_count:
                break
                
            summary_handle = Entrez.esummary(db="nucleotide", id=ncb_id)
            summary_record = Entrez.read(summary_handle)
            summary_handle.close()
            
            # Extract out the dictionary safely even if NCBI returns it inside a ListElement container
            if isinstance(summary_record, list) or hasattr(summary_record, '__getitem__') and not hasattr(summary_record, 'get'):
                record_data = summary_record[0]
            else:
                record_data = summary_record
                
            accession = record_data.get("Caption", ncb_id)
            title = record_data.get("Title", f"Phage_{ncb_id}")
            
            if accession in seen_accessions or "_" in accession:
                continue
                
            seen_accessions.add(accession)
            print(f"📥 [{downloaded_count + 1}/{target_count}] Pulling sequence {accession} ({title[:40]}...)")
            
            fetch_handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            fasta_data = fetch_handle.read()
            fetch_handle.close()
            
            nested_folder = os.path.join(output_dir, accession)
            os.makedirs(nested_folder, exist_ok=True)
            
            file_path = os.path.join(nested_folder, f"{accession}.fasta")
            with open(file_path, "w") as f:
                f.write(fasta_data)
                
            downloaded_count += 1
            time.sleep(0.5)
            
        print(f"\n✅ Successfully downloaded {downloaded_count} unique Klebsiella phage genomes!")

    except Exception as e:
        print(f"❌ An error occurred during retrieval: {e}")

if __name__ == "__main__":
    fetch_unique_klebsiella_phages(target_count=10)
