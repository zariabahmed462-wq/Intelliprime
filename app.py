"""
IntelliPrime - Web Interface
Powered by Streamlit
"""

import streamlit as st
from engine import fetch_gene, extract_cds, codon_optimize, check_internal_sites, add_cloning_sites, run_primer_design

st.set_page_config(page_title="IntelliPrime", page_icon="🧬")
st.title("🧬 IntelliPrime")
st.subheader("Context-Aware Primer Design Engine")

# Simple interface
accession = st.text_input("Enter Accession Number:", value="NM_000518.5")

goal = st.selectbox("Select Goal:", ["Cloning into pET28a (NdeI/XhoI)"])

if st.button("🚀 Design Primers"):
    with st.spinner("Fetching gene from NCBI..."):
        try:
            record = fetch_gene(accession)
            cds = extract_cds(record)
            
            st.success(f"Gene loaded: {len(cds)} bp CDS")
            
            with st.expander("View original CDS (first 100bp)"):
                st.code(cds[:100] + "...")
            
            # Codon optimize
            optimized = codon_optimize(cds)
            
            with st.expander("View optimized CDS (first 100bp)"):
                st.code(optimized[:100] + "...")
            
            # Check internal sites
            problems = check_internal_sites(optimized)
            if problems:
                for p in problems:
                    st.warning(p)
            else:
                st.success("No internal restriction sites found")
            
            # Add cloning sites
            template = add_cloning_sites(optimized)
            
            # Design primers
            target_start = len("GATATACATATG")
            primers = run_primer_design(template, target_start, len(optimized))
            
            if primers:
                st.subheader("🎯 Designed Primers")
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Forward Primer", primers['forward'])
                    st.caption(f"Tm: {primers['fwd_tm']:.1f}°C")
                with col2:
                    st.metric("Reverse Primer", primers['reverse'])
                    st.caption(f"Tm: {primers['rev_tm']:.1f}°C")
                st.info(f"Product Size: {primers['product_size']} bp")
            
            st.subheader("📋 Complete Template")
            st.code(template, language="text")
            
        except Exception as e:
            st.error(f"Error: {str(e)}")

st.caption("IntelliPrime v0.1 | Open Source | MIT License")