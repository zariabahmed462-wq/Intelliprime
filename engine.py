"""
IntelliPrime - Context-Aware Primer Design Engine
Week 1: Core Engine - Cloning Goal (pET28a NdeI/XhoI)
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import primer3
from python_codon_tables import get_codons_table

# Always identify yourself to NCBI
Entrez.email = "your_email@example.com"


def fetch_gene(accession):
    """
    Fetch a GenBank record by accession number.
    Returns the SeqRecord object.
    """
    print(f"Fetching {accession} from NCBI...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    print(f"Loaded: {record.annotations.get('organism', 'Unknown')} - {len(record.seq)} bp")
    return record


def extract_cds(record):
    """
    Extract the CDS (coding sequence) feature from the GenBank record.
    Returns the CDS sequence as a string.
    """
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            print(f"CDS found: {len(cds_seq)} bp")
            return str(cds_seq)
    
    print("Warning: No CDS feature found. Using full sequence.")
    return str(record.seq)


def codon_optimize(dna_sequence):
    """
    Optimize DNA sequence for E. coli K12 using most frequent codons.
    Returns optimized DNA string.
    """
    codon_table = get_codons_table("e_coli")
    
    # Translate DNA to protein
    seq_obj = Seq(dna_sequence)
    protein = str(seq_obj.translate())
    
    # Reverse translate using most frequent codon for each amino acid
    optimized = ""
    for aa in protein:
        if aa == "*":  # Stop codon
            optimized += "TAA"  # Most common E. coli stop codon
        else:
            # Get the codon with highest frequency
            best_codon = max(codon_table[aa], key=codon_table[aa].get)
            optimized += best_codon
    
    print(f"Codon optimized: {len(optimized)} bp")
    return optimized


def check_internal_sites(sequence):
    """
    Check for internal restriction sites that would interfere with cloning.
    Returns list of problems found.
    """
    problems = []
    
    # Restriction sites we're using
    sites_to_check = {
        "NdeI": "CATATG",
        "XhoI": "CTCGAG",
        "NcoI": "CCATGG",  # Common in pET vectors
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
    }
    
    for name, site in sites_to_check.items():
        if site in sequence:
            problems.append(f"WARNING: Internal {name} site ({site}) found!")
    
    return problems


def add_cloning_sites(optimized_cds):
    """
    Add 5' clamp + NdeI and 3' XhoI + clamp for pET28a cloning.
    Returns the complete template sequence.
    """
    five_prime = "GATATA"   # Random clamp for efficient cutting
    ndei_site = "CATATG"    # NdeI - contains start codon ATG
    three_prime = "CTCGAG"  # XhoI
    clamp_3 = "GTG"         # Random 3' clamp
    
    complete_template = five_prime + ndei_site + optimized_cds + three_prime + clamp_3
    
    print(f"Final template: {len(complete_template)} bp")
    print(f"  5' clamp + NdeI: {five_prime}{ndei_site}")
    print(f"  3' XhoI + clamp: {three_prime}{clamp_3}")
    
    return complete_template


def run_primer_design(template, target_start, target_length):
    """
    Run Primer3 to design primers for amplifying the target region.
    Returns primer pair dictionary or None.
    """
    # The region we want to amplify (the coding sequence between the RE sites)
    seq_args = {
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_INCLUDED_REGION": [target_start, target_length],
    }
    
    # Design parameters for cloning
    global_args = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 58.0,
        "PRIMER_MAX_TM": 62.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[target_length - 20, target_length + 20]],
    }
    
    # Run Primer3
    result = primer3.bindings.design_primers(seq_args, global_args)
    
    if result["PRIMER_PAIR_NUM_RETURNED"] > 0:
        return {
            "forward": result["PRIMER_LEFT_0_SEQUENCE"],
            "reverse": result["PRIMER_RIGHT_0_SEQUENCE"],
            "fwd_tm": result["PRIMER_LEFT_0_TM"],
            "rev_tm": result["PRIMER_RIGHT_0_TM"],
            "product_size": result["PRIMER_PAIR_0_PRODUCT_SIZE"],
        }
    
    return None


def main(accession):
    """
    Main pipeline: Fetch gene → Extract CDS → Optimize → Add cloning sites → Design primers.
    """
    print("\n" + "="*60)
    print(f"IntelliPrime - Goal: Cloning into pET28a (NdeI/XhoI)")
    print(f"Target Accession: {accession}")
    print("="*60 + "\n")
    
    # Step 1: Fetch the gene
    record = fetch_gene(accession)
    
    # Step 2: Extract CDS
    cds = extract_cds(record)
    print(f"Original CDS (first 60bp): {cds[:60]}...\n")
    
    # Step 3: Codon optimize
    optimized = codon_optimize(cds)
    print(f"Optimized (first 60bp): {optimized[:60]}...\n")
    
    # Step 4: Check for internal sites in the optimized sequence
    print("Checking for internal restriction sites...")
    problems = check_internal_sites(optimized)
    for problem in problems:
        print(f"  {problem}")
    if not problems:
        print("  No conflicts found.\n")
    
    # Step 5: Add cloning sites
    template = add_cloning_sites(optimized)
    print(f"Complete template (first 80bp): {template[:80]}...")
    print(f"Complete template (last 80bp): ...{template[-80:]}\n")
    
    # Step 6: Design primers
    print("Designing primers with Primer3...")
    # Target region starts after 5' clamp+NdeI, length is the optimized CDS
    target_start = len("GATATACATATG")
    target_length = len(optimized)
    
    primers = run_primer_design(template, target_start, target_length)
    
    if primers:
        print("Primers designed successfully:")
        print(f"  Forward ({primers['fwd_tm']:.1f}°C): {primers['forward']}")
        print(f"  Reverse ({primers['rev_tm']:.1f}°C): {primers['reverse']}")
        print(f"  Product size: {primers['product_size']} bp")
    else:
        print("Primer3 could not find suitable primers with current constraints.")
    
    print("\n" + "="*60)
    print("Design complete!")
    print("="*60)
    
    return {
        "accession": accession,
        "cds_length": len(cds),
        "optimized_cds": optimized,
        "complete_template": template,
        "primers": primers,
        "warnings": problems,
    }


if __name__ == "__main__":
    # Test with Human Hemoglobin Beta (HBB)
    result = main("NM_000518.5")