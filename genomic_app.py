# app/genomic_app.py
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import tempfile
import os

# Import backend
from app.backend import sequence_analyzer, variant_analyzer, genomic_utilities

# --- Streamlit Config ---
st.set_page_config(
    page_title="Genomic Analysis Toolkit",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- Custom CSS ---
st.markdown("""
<style>
    .main-header { font-size: 2.5rem; color: #2E86AB; text-align: center; margin-bottom: 2rem; font-weight: bold; }
    .section-header { font-size: 1.8rem; color: #A23B72; margin-top: 2rem; margin-bottom: 1rem; border-bottom: 2px solid #A23B72; padding-bottom: 0.5rem; }
</style>
""", unsafe_allow_html=True)

# --- Main App ---
def main():
    st.markdown('<h1 class="main-header">🧬 Genomic Analysis Toolkit</h1>', unsafe_allow_html=True)

    with st.sidebar:
        st.markdown("## 🔍 Navigation")
        app_mode = st.radio(
            "Choose Analysis Type:",
            ["🏠 Home", "🧬 Sequence Analysis", "📊 Variant Analysis", "ℹ️ About"]
        )

    if "Home" in app_mode:
        show_home()
    elif "Sequence Analysis" in app_mode:
        sequence_analysis()
    elif "Variant Analysis" in app_mode:
        variant_analysis()
    elif "About" in app_mode:
        show_about()


# --- Home ---
def show_home():
    st.markdown("## Welcome to the Genomic Analysis Toolkit!")
    st.markdown("""
    This application provides:
    - **🧬 Sequence Analysis** (FASTA)
    - **📊 Variant Analysis** (VCF)
    - **Utility tools** like GC calculation, MW, melting temp, etc.
    """)


# --- Sequence Analysis ---
def sequence_analysis():
    st.markdown('<h2 class="section-header">🧬 Sequence Analysis</h2>', unsafe_allow_html=True)

    uploaded_file = st.file_uploader("📤 Upload FASTA File", type=['fasta', 'fa', 'fna', 'txt'])

    if uploaded_file is not None:
        success, msg = sequence_analyzer.load_sequences(uploaded_file, "fasta")

        if success:
            st.success(msg)

            seq_ids = list(sequence_analyzer.sequences.keys())
            selected_seq_id = st.selectbox("Select sequence for analysis:", seq_ids)
            sequence = sequence_analyzer.sequences[selected_seq_id]

            # --- Stats ---
            stats = sequence_analyzer.get_sequence_stats(sequence)
            df_stats = pd.DataFrame([stats])
            st.dataframe(df_stats)

            # --- GC Sliding Plot ---
            window = st.slider("Window size for GC calculation", 50, 1000, 100, 50)
            if st.button("Generate GC Plot"):
                positions, gc_values = sequence_analyzer.calculate_gc_content_sliding(sequence, window)
                fig, ax = plt.subplots(figsize=(10, 5))
                ax.plot(positions, gc_values, color='#2E86AB', linewidth=2)
                ax.set_xlabel("Position (bp)")
                ax.set_ylabel("GC Content (%)")
                ax.set_title(f"GC Content Sliding Window - {selected_seq_id}")
                st.pyplot(fig)

            # --- Motif Search ---
            motif = st.text_input("Find motif (e.g. ATG):")
            if motif:
                positions = sequence_analyzer.find_motifs(sequence, motif)
                st.info(f"Found {len(positions)} occurrences of '{motif}'")
                st.write(positions)
        else:
            st.error(msg)
    else:
        st.warning("👆 Upload a FASTA file to begin.")


# --- Variant Analysis ---
def variant_analysis():
    st.markdown('<h2 class="section-header">📊 Variant Analysis</h2>', unsafe_allow_html=True)
    uploaded_file = st.file_uploader("📤 Upload VCF File", type=['vcf', 'vcf.gz'])

    if uploaded_file is not None:
        # Save temporarily to disk for cyvcf2
        with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as tmp:
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name

        success, msg = variant_analyzer.load_vcf_file(tmp_path)
        os.remove(tmp_path)

        if success:
            st.success(msg)
            summary = variant_analyzer.get_variant_summary()
            st.dataframe(pd.DataFrame([summary]))

            # --- Chromosome Distribution ---
            chrom_counts = variant_analyzer.get_chromosome_distribution()
            fig, ax = plt.subplots(figsize=(8, 5))
            chrom_counts.plot(kind='bar', color='#FF6B6B', ax=ax)
            ax.set_title("Variants per Chromosome")
            ax.set_ylabel("Count")
            st.pyplot(fig)
        else:
            st.error(msg)
    else:
        st.warning("👆 Upload a VCF file to analyze.")


# --- About ---
def show_about():
    st.markdown('<h2 class="section-header">ℹ️ About</h2>', unsafe_allow_html=True)
    st.markdown("""
    **Genomic Analysis Toolkit** combines Biopython, Pandas, and Streamlit
    to enable interactive genomic data exploration.
    """)


# --- Run App ---
if __name__ == "__main__":
    main()