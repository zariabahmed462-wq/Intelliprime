"""
IntelliPrime - Web Interface
Powered by Streamlit
"""

import streamlit as st
from engine import (
    fetch_gene, extract_cds, codon_optimize,
    check_internal_sites, add_cloning_sites
)

st.set_page_config(page_title="IntelliPrime", page_icon="🧬")
st.title("🧬 IntelliPrime")
st.subheader("Context-Aware Primer Design Engine")

accession = st.text_input("Enter Accession Number:", value="NM_000518.5")
goal = st.selectbox("Select Goal:", ["Cloning into pET28a (NdeI/XhoI)"])

if st.button("🚀 Design Primers"):
    with st.spinner("Fetching gene from NCBI..."):
        try:
            record = fetch_gene(accession)
            cds = extract_cds(record)
            
            st.success(f"Gene loaded: {len(cds)} bp CDS")
            st.write(f"**Organism:** {record.annotations.get('organism', 'Unknown')}")
            
            with st.expander("🔍 View original CDS (first 100bp)"):
                st.code(cds[:100] + "...")
            
            optimized = codon_optimize(cds)
            
            with st.expander("🔍 View codon-optimized CDS (first 100bp)"):
                st.code(optimized[:100] + "...")
            
            problems = check_internal_sites(optimized)
            if problems:
                for p in problems:
                    st.warning(p)
            else:
                st.success("✅ No internal restriction sites found")
            
            template = add_cloning_sites(optimized)
            
            st.subheader("🎯 Designed Primers (Simple Mode)")
            st.info("Using simple design rules. Primer3 coming soon.")
            
            # Simple primer design
            fwd_primer = template[:25]
            rev_primer = template[-25:]
            
            col1, col2 = st.columns(2)
            with col1:
                st.text_area("Forward Primer", fwd_primer, height=80)
            with col2:
                st.text_area("Reverse Primer", rev_primer, height=80)
            
            st.subheader("📋 Complete Template Sequence")
            st.text_area("Template", template, height=100)
            st.caption(f"Total length: {len(template)} bp")
            
            st.download_button(
                "📥 Download Template (FASTA)",
                f">IntelliPrime_Cloning_Template_{accession}\n{template}",
                file_name=f"{accession}_template.fasta"
            )
            
        except Exception as e:
            st.error(f"Error: {str(e)}")

st.divider()
st.caption("🔬 IntelliPrime v0.1 | Open Source | MIT License")
st.caption("Developed as a research tool for intelligent primer design")