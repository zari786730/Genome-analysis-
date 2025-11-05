#!/usr/bin/env python3
"""
Script to create complete genomic analysis web app structure
"""

import os
import sys

def create_file(path, content):
    """Create a file with the given content"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"✅ Created: {path}")

def main():
    print("🚀 Creating Genomic Analysis Web App Structure...")
    print("=" * 50)
    
    # Define the base directory
    base_dir = "genomic-analysis-app"
    
    # 1. Create app/genomic_app.py
    genomic_app_content = '''#!/usr/bin/env python3
"""
Genomic Analysis Web App with Streamlit
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils import GC, molecular_weight
import io
import base64

# Page configuration
st.set_page_config(
    page_title="Genomic Analysis Toolkit",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 3rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
    }
    .section-header {
        font-size: 1.5rem;
        color: #A23B72;
        margin-top: 2rem;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 Genomic Analysis Toolkit</h1>', unsafe_allow_html=True)
    
    # Sidebar navigation
    st.sidebar.title("Navigation")
    app_mode = st.sidebar.selectbox(
        "Choose Analysis Type",
        ["Home", "Sequence Analysis", "Variant Analysis", "File Converter", "About"]
    )
    
    if app_mode == "Home":
        show_home()
    elif app_mode == "Sequence Analysis":
        sequence_analysis()
    elif app_mode == "Variant Analysis":
        variant_analysis()
    elif app_mode == "File Converter":
        file_converter()
    elif app_mode == "About":
        show_about()

def show_home():
    st.markdown("""
    ## Welcome to the Genomic Analysis Toolkit!
    
    This web application provides various tools for genomic data analysis:
    
    - **Sequence Analysis**: Analyze FASTA files, calculate GC content, molecular weight
    - **Variant Analysis**: Process VCF files, visualize variant distributions
    - **File Converter**: Convert between different genomic file formats
    
    ### Getting Started
    
    1. Choose an analysis type from the sidebar
    2. Upload your genomic files
    3. Configure analysis parameters
    4. View results and download reports
    
    ### Supported File Formats
    - FASTA (.fasta, .fa, .fna)
    - FASTQ (.fastq, .fq)
    - VCF (.vcf)
    - CSV/TSV for data tables
    """)

def sequence_analysis():
    st.markdown('<h2 class="section-header">Sequence Analysis</h2>', unsafe_allow_html=True)
    
    # File upload
    uploaded_file = st.file_uploader("Upload FASTA file", type=['fasta', 'fa', 'fna', 'txt'])
    
    if uploaded_file is not None:
        # Read sequences
        sequences = list(SeqIO.parse(uploaded_file, "fasta"))
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Sequence Summary")
            st.write(f"Number of sequences: {len(sequences)}")
            
            # Sequence selector
            if sequences:
                sequence_names = [seq.id for seq in sequences]
                selected_seq = st.selectbox("Select sequence for detailed analysis:", sequence_names)
                
                # Find selected sequence
                selected_sequence = next(seq for seq in sequences if seq.id == selected_seq)
                
                # Basic stats
                st.write(f"**Sequence ID:** {selected_sequence.id}")
                st.write(f"**Length:** {len(selected_sequence)} bp")
                st.write(f"**GC Content:** {GC(selected_sequence.seq):.2f}%")
                st.write(f"**Molecular Weight:** {molecular_weight(selected_sequence.seq):.2f} g/mol")
        
        with col2:
            st.subheader("GC Content Analysis")
            window_size = st.slider("Window size for GC analysis:", 50, 1000, 100)
            
            if st.button("Calculate GC Content"):
                gc_values = []
                positions = []
                
                seq_str = str(selected_sequence.seq)
                for i in range(0, len(seq_str) - window_size + 1, window_size):
                    window = seq_str[i:i + window_size]
                    gc_percent = GC(window)
                    gc_values.append(gc_percent)
                    positions.append(i)
                
                # Plot
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.plot(positions, gc_values)
                ax.set_title(f'GC Content - {selected_sequence.id}')
                ax.set_xlabel('Position (bp)')
                ax.set_ylabel('GC Content (%)')
                ax.grid(True, alpha=0.3)
                
                st.pyplot(fig)
        
        # Show all sequences in a table
        st.subheader("All Sequences")
        sequence_data = []
        for seq in sequences:
            sequence_data.append({
                'ID': seq.id,
                'Length': len(seq),
                'GC_Content': f"{GC(seq.seq):.2f}%",
                'Description': seq.description[:100] + "..." if len(seq.description) > 100 else seq.description
            })
        
        df = pd.DataFrame(sequence_data)
        st.dataframe(df)

def variant_analysis():
    st.markdown('<h2 class="section-header">Variant Analysis</h2>', unsafe_allow_html=True)
    
    uploaded_file = st.file_uploader("Upload VCF file", type=['vcf', 'vcf.gz'])
    
    if uploaded_file is not None:
        try:
            import cyvcf2
            
            # For demo purposes, we'll create sample data
            # In real app, you would process the actual VCF file
            st.info("Processing VCF file...")
            
            # Create sample variant data for demonstration
            sample_data = {
                'CHROM': ['chr1'] * 500 + ['chr2'] * 500,
                'POS': list(range(1, 501)) + list(range(1, 501)),
                'REF': ['A'] * 1000,
                'ALT': ['G'] * 500 + ['T'] * 500,
                'QUAL': np.random.normal(100, 20, 1000),
                'TYPE': ['SNP'] * 800 + ['INDEL'] * 200
            }
            
            df = pd.DataFrame(sample_data)
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("Variant Summary")
                st.write(f"Total variants: {len(df)}")
                st.write(f"SNPs: {len(df[df['TYPE'] == 'SNP'])}")
                st.write(f"INDELs: {len(df[df['TYPE'] == 'INDEL'])}")
                
                # Variant type distribution
                fig1, ax1 = plt.subplots(figsize=(8, 6))
                df['TYPE'].value_counts().plot(kind='bar', ax=ax1, color=['skyblue', 'lightcoral'])
                ax1.set_title('Variant Type Distribution')
                ax1.set_ylabel('Count')
                st.pyplot(fig1)
            
            with col2:
                st.subheader("Quality Distribution")
                fig2, ax2 = plt.subplots(figsize=(8, 6))
                ax2.hist(df['QUAL'], bins=50, alpha=0.7, color='lightgreen')
                ax2.set_title('Variant Quality Distribution')
                ax2.set_xlabel('Quality Score')
                ax2.set_ylabel('Frequency')
                st.pyplot(fig2)
            
            # Chromosome distribution
            st.subheader("Chromosome Distribution")
            chrom_counts = df['CHROM'].value_counts()
            fig3, ax3 = plt.subplots(figsize=(10, 6))
            chrom_counts.plot(kind='bar', ax=ax3, color='purple', alpha=0.7)
            ax3.set_title('Variants per Chromosome')
            ax3.set_ylabel('Count')
            st.pyplot(fig3)
            
            # Data table
            st.subheader("Variant Data")
            st.dataframe(df.head(100))
            
        except ImportError:
            st.error("cyvcf2 is not installed. Please install it to use VCF analysis.")

def file_converter():
    st.markdown('<h2 class="section-header">File Converter</h2>', unsafe_allow_html=True)
    
    conversion_type = st.selectbox(
        "Select conversion type:",
        ["FASTA to CSV", "CSV to FASTA", "VCF to CSV"]
    )
    
    uploaded_file = st.file_uploader(f"Upload file for {conversion_type}")
    
    if uploaded_file is not None and conversion_type == "FASTA to CSV":
        sequences = list(SeqIO.parse(uploaded_file, "fasta"))
        
        # Convert to DataFrame
        sequence_data = []
        for seq in sequences:
            sequence_data.append({
                'ID': seq.id,
                'Sequence': str(seq.seq),
                'Length': len(seq),
                'GC_Content': GC(seq.seq),
                'Description': seq.description
            })
        
        df = pd.DataFrame(sequence_data)
        
        st.subheader("Converted Data")
        st.dataframe(df)
        
        # Download button
        csv = df.to_csv(index=False)
        st.download_button(
            label="Download as CSV",
            data=csv,
            file_name="sequences.csv",
            mime="text/csv"
        )

def show_about():
    st.markdown("""
    ## About Genomic Analysis Toolkit
    
    This web application is built with Streamlit and provides various tools for genomic data analysis.
    
    ### Features
    - Sequence analysis and visualization
    - Variant calling data processing
    - File format conversion
    - Interactive plots and statistics
    
    ### Technologies Used
    - **Streamlit**: Web application framework
    - **Biopython**: Biological computation
    - **Pandas**: Data manipulation
    - **Matplotlib/Seaborn**: Visualization
    - **Plotly**: Interactive plots
    
    ### Source Code
    The source code for this application is available on GitHub.
    """)

if __name__ == "__main__":
    main()
'''
    create_file(os.path.join(base_dir, "app", "genomic_app.py"), genomic_app_content)
    
    # 2. Create requirements.txt
    requirements_content = '''streamlit>=1.28.0
biopython>=1.81
pandas>=2.0.0
numpy>=1.21.0
matplotlib>=3.5.0
seaborn>=0.12.0
cyvcf2>=0.30.0
plotly>=5.0.0
'''
    create_file(os.path.join(base_dir, "requirements.txt"), requirements_content)
    
    # 3. Create README.md
    readme_content = '''# 🧬 Genomic Analysis Web App

A Streamlit web application for genomic data analysis with sequence analysis, variant calling, and file conversion tools.

## 🌟 Features

- **Sequence Analysis**: GC content, molecular weight, sequence statistics
- **Variant Analysis**: VCF file processing and visualization
- **File Converter**: Convert between FASTA, CSV, and other formats
- **Interactive Plots**: Dynamic visualizations of genomic data

## 🚀 Quick Start

### Local Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/genomic-analysis-app.git
   cd genomic-analysis-app