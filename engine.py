"""
IntelliPrime - Context-Aware Primer Design Engine
Week 1: Core Engine - Cloning Goal (pET28a NdeI/XhoI)
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from python_codon_tables import get_codons_table

Entrez.email = "your_email@example.com"


def fetch_gene(accession):
    """Fetch a GenBank record by accession number."""
    print(f"Fetching {accession} from NCBI...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    print(f"Loaded: {record.annotations.get('organism', 'Unknown')} - {len(record.seq)} bp")
    return record


def extract_cds(record):
    """Extract the CDS feature from the GenBank record."""
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            print(f"CDS found: {len(cds_seq)} bp")
            return str(cds_seq)
    
    print("Warning: No CDS feature found. Using full sequence.")
    return str(record.seq)


def codon_optimize(dna_sequence):
    """Optimize DNA sequence for E. coli K12 using most frequent codons."""
    codon_table = get_codons_table("e_coli")
    
    seq_obj = Seq(dna_sequence)
    protein = str(seq_obj.translate())
    
    optimized = ""
    for aa in protein:
        if aa == "*":
            optimized += "TAA"
        else:
            best_codon = max(codon_table[aa], key=codon_table[aa].get)
            optimized += best_codon
    
    print(f"Codon optimized: {len(optimized)} bp")
    return optimized


def check_internal_sites(sequence):
    """Check for internal restriction sites."""
    problems = []
    
    sites_to_check = {
        "NdeI": "CATATG",
        "XhoI": "CTCGAG",
        "NcoI": "CCATGG",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
    }
    
    for name, site in sites_to_check.items():
        if site in sequence:
            problems.append(f"WARNING: Internal {name} site ({site}) found!")
    
    return problems


def add_cloning_sites(optimized_cds):
    """Add 5' clamp + NdeI and 3' XhoI + clamp for pET28a cloning."""
    five_prime = "GATATA"
    ndei_site = "CATATG"
    three_prime = "CTCGAG"
    clamp_3 = "GTG"
    
    complete_template = five_prime + ndei_site + optimized_cds + three_prime + clamp_3
    
    print(f"Final template: {len(complete_template)} bp")
    return complete_template


def main(accession):
    """Main pipeline."""
    print("\n" + "="*60)
    print(f"IntelliPrime - Goal: Cloning into pET28a (NdeI/XhoI)")
    print(f"Target Accession: {accession}")
    print("="*60 + "\n")
    
    record = fetch_gene(accession)
    cds = extract_cds(record)
    optimized = codon_optimize(cds)
    
    problems = check_internal_sites(optimized)
    for problem in problems:
        print(f"  {problem}")
    
    template = add_cloning_sites(optimized)
    
    print("\n" + "="*60)
    print("Design complete!")
    print("="*60)
    
    return {
        "accession": accession,
        "cds_length": len(cds),
        "optimized_cds": optimized,
        "complete_template": template,
        "warnings": problems,
    }


if __name__ == "__main__":
    result = main("NM_000518.5")