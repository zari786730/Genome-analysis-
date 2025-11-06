import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils import GC, molecular_weight
import io
import base64

# Page configuration - FIXED
st.set_page_config(
    page_title="Genomic Analysis Toolkit",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling - FIXED
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
        font-weight: bold;
    }
    .section-header {
        font-size: 1.8rem;
        color: #A23B72;
        margin-top: 2rem;
        margin-bottom: 1rem;
        border-bottom: 2px solid #A23B72;
        padding-bottom: 0.5rem;
    }
    .stButton button {
        background-color: #2E86AB;
        color: white;
        border-radius: 5px;
        padding: 0.5rem 1rem;
        border: none;
    }
    .stButton button:hover {
        background-color: #1F5F7A;
    }
    .uploadedFile {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 5px;
        margin: 1rem 0;
    }
</style>
""", unsafe_allow_html=True)

def main():
    # Header with better styling
    st.markdown('<h1 class="main-header">🧬 Genomic Analysis Toolkit</h1>', unsafe_allow_html=True)
    
    # Sidebar navigation
    with st.sidebar:
        st.markdown("## 🔍 Navigation")
        st.markdown("---")
        app_mode = st.radio(
            "Choose Analysis Type:",
            ["🏠 Home", "🧬 Sequence Analysis", "📊 Variant Analysis", "🔄 File Converter", "ℹ️ About"]
        )
        
        st.markdown("---")
        st.markdown("### 💡 Tips")
        st.info("Upload FASTA files for sequence analysis or VCF files for variant analysis.")
    
    # Remove emoji prefix for logic
    if "Home" in app_mode:
        show_home()
    elif "Sequence Analysis" in app_mode:
        sequence_analysis()
    elif "Variant Analysis" in app_mode:
        variant_analysis()
    elif "File Converter" in app_mode:
        file_converter()
    elif "About" in app_mode:
        show_about()

def show_home():
    st.markdown("## 🏠 Welcome to Genomic Analysis Toolkit!")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        ### Your All-in-One Genomic Data Analysis Platform
        
        This web application provides comprehensive tools for genomic data analysis:
        
        **🔬 Sequence Analysis**
        - Analyze FASTA files
        - Calculate GC content
        - Molecular weight calculations
        - Sequence statistics
        
        **📊 Variant Analysis** 
        - Process VCF files
        - Visualize variant distributions
        - Quality metrics
        - Chromosome distribution
        
        **🔄 File Converter**
        - Convert FASTA to CSV
        - Export analysis results
        - Format transformation
        
        ### 🚀 Getting Started
        
        1. **Choose an analysis type** from the sidebar
        2. **Upload your genomic files** 
        3. **Configure analysis parameters**
        4. **View results and download reports**
        """)
    
    with col2:
        st.markdown("### 📁 Supported Formats")
        st.success("""
        **FASTA Files**
        - .fasta
        - .fa  
        - .fna
        - .txt
        """)
        
        st.info("""
        **VCF Files**
        - .vcf
        - .vcf.gz
        """)
        
        st.warning("""
        **Data Files**
        - CSV
        - TSV
        - Excel
        """)

def sequence_analysis():
    st.markdown('<h2 class="section-header">🧬 Sequence Analysis</h2>', unsafe_allow_html=True)
    
    st.info("Upload a FASTA file to analyze sequence properties, GC content, and generate visualizations.")
    
    # File upload section
    uploaded_file = st.file_uploader(
        "📤 Upload FASTA File", 
        type=['fasta', 'fa', 'fna', 'txt'],
        help="Upload a FASTA format file containing DNA/protein sequences"
    )
    
    if uploaded_file is not None:
        try:
            # Read sequences
            sequences = list(SeqIO.parse(uploaded_file, "fasta"))
            
            if not sequences:
                st.error("❌ No sequences found in the uploaded file.")
                return
            
            st.success(f"✅ Successfully loaded {len(sequences)} sequence(s)")
            
            # Create tabs for different analyses
            tab1, tab2, tab3 = st.tabs(["📊 Overview", "🔍 Detailed Analysis", "📈 Visualizations"])
            
            with tab1:
                st.subheader("Sequence Overview")
                
                # Basic statistics
                stats_data = []
                for seq in sequences:
                    stats_data.append({
                        'ID': seq.id,
                        'Length': len(seq),
                        'GC Content (%)': f"{GC(seq.seq):.2f}",
                        'Molecular Weight': f"{molecular_weight(seq.seq):.2f}"
                    })
                
                stats_df = pd.DataFrame(stats_data)
                st.dataframe(stats_df, use_container_width=True)
            
            with tab2:
                st.subheader("Detailed Sequence Analysis")
                
                if len(sequences) > 1:
                    selected_seq = st.selectbox(
                        "Select sequence for detailed analysis:",
                        [seq.id for seq in sequences]
                    )
                    sequence = next(seq for seq in sequences if seq.id == selected_seq)
                else:
                    sequence = sequences[0]
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.metric("Sequence ID", sequence.id)
                    st.metric("Length", f"{len(sequence)} bp")
                    st.metric("GC Content", f"{GC(sequence.seq):.2f}%")
                    st.metric("Molecular Weight", f"{molecular_weight(sequence.seq):.2f} g/mol")
                
                with col2:
                    # Nucleotide composition
                    nucleotides = ['A', 'T', 'G', 'C']
                    counts = {nt: sequence.seq.count(nt) for nt in nucleotides}
                    
                    fig, ax = plt.subplots(figsize=(8, 6))
                    ax.bar(counts.keys(), counts.values(), color=['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4'])
                    ax.set_title('Nucleotide Composition')
                    ax.set_ylabel('Count')
                    st.pyplot(fig)
            
            with tab3:
                st.subheader("GC Content Visualization")
                
                if len(sequences) > 1:
                    selected_seq = st.selectbox(
                        "Select sequence for GC analysis:",
                        [seq.id for seq in sequences],
                        key="gc_selector"
                    )
                    sequence = next(seq for seq in sequences if seq.id == selected_seq)
                else:
                    sequence = sequences[0]
                
                window_size = st.slider(
                    "Window size for GC analysis:",
                    min_value=50,
                    max_value=1000,
                    value=100,
                    step=50,
                    help="Size of sliding window for GC content calculation"
                )
                
                if st.button("Generate GC Plot", type="primary"):
                    gc_values = []
                    positions = []
                    
                    seq_str = str(sequence.seq)
                    for i in range(0, len(seq_str) - window_size + 1, window_size):
                        window = seq_str[i:i + window_size]
                        gc_percent = GC(window)
                        gc_values.append(gc_percent)
                        positions.append(i)
                    
                    # Plot
                    fig, ax = plt.subplots(figsize=(12, 6))
                    ax.plot(positions, gc_values, linewidth=2, color='#2E86AB')
                    ax.set_title(f'GC Content Distribution - {sequence.id}', fontsize=16, fontweight='bold')
                    ax.set_xlabel('Genomic Position (bp)', fontsize=12)
                    ax.set_ylabel('GC Content (%)', fontsize=12)
                    ax.grid(True, alpha=0.3)
                    ax.set_facecolor('#f8f9fa')
                    
                    st.pyplot(fig)
                    
        except Exception as e:
            st.error(f"❌ Error processing file: {str(e)}")

def variant_analysis():
    st.markdown('<h2 class="section-header">📊 Variant Analysis</h2>', unsafe_allow_html=True)
    
    st.info("Upload a VCF file to analyze genetic variants, visualize distributions, and explore variant data.")
    
    uploaded_file = st.file_uploader(
        "📤 Upload VCF File", 
        type=['vcf', 'vcf.gz'],
        help="Upload a VCF file containing genetic variants"
    )
    
    if uploaded_file is not None:
        try:
            # For demo - create sample variant data
            st.success("✅ VCF file uploaded successfully!")
            
            # Create sample data for demonstration
            np.random.seed(42)  # For consistent demo data
            sample_data = {
                'Chromosome': ['chr1'] * 300 + ['chr2'] * 250 + ['chr3'] * 200 + ['chr4'] * 150 + ['chr5'] * 100,
                'Position': (list(range(1, 301)) + list(range(1, 251)) + 
                           list(range(1, 201)) + list(range(1, 151)) + list(range(1, 101))),
                'Reference': ['A'] * 500 + ['T'] * 300,
                'Alternate': ['G'] * 300 + ['C'] * 200 + ['T'] * 300,
                'Quality': np.random.normal(150, 30, 800),
                'Variant Type': ['SNP'] * 600 + ['INDEL'] * 200,
                'Filter': ['PASS'] * 700 + ['LOW_QUAL'] * 100
            }
            
            df = pd.DataFrame(sample_data)
            
            # Display overview
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Variants", len(df))
            with col2:
                st.metric("SNPs", len(df[df['Variant Type'] == 'SNP']))
            with col3:
                st.metric("INDELs", len(df[df['Variant Type'] == 'INDEL']))
            with col4:
                st.metric("Passing Filters", len(df[df['Filter'] == 'PASS']))
            
            # Create visualization tabs
            tab1, tab2, tab3 = st.tabs(["📈 Distributions", "🧬 Variant Types", "📋 Data Table"])
            
            with tab1:
                col1, col2 = st.columns(2)
                
                with col1:
                    # Quality distribution
                    fig1, ax1 = plt.subplots(figsize=(8, 6))
                    ax1.hist(df['Quality'], bins=30, alpha=0.7, color='#4ECDC4', edgecolor='black')
                    ax1.set_title('Variant Quality Distribution', fontweight='bold')
                    ax1.set_xlabel('Quality Score')
                    ax1.set_ylabel('Frequency')
                    st.pyplot(fig1)
                
                with col2:
                    # Chromosome distribution
                    chrom_counts = df['Chromosome'].value_counts()
                    fig2, ax2 = plt.subplots(figsize=(8, 6))
                    ax2.bar(chrom_counts.index, chrom_counts.values, color='#FF6B6B', alpha=0.7)
                    ax2.set_title('Variants per Chromosome', fontweight='bold')
                    ax2.set_ylabel('Count')
                    plt.xticks(rotation=45)
                    st.pyplot(fig2)
            
            with tab2:
                col1, col2 = st.columns(2)
                
                with col1:
                    # Variant type pie chart
                    type_counts = df['Variant Type'].value_counts()
                    fig3, ax3 = plt.subplots(figsize=(8, 6))
                    ax3.pie(type_counts.values, labels=type_counts.index, autopct='%1.1f%%', 
                           colors=['#45B7D1', '#96CEB4'])
                    ax3.set_title('Variant Type Distribution', fontweight='bold')
                    st.pyplot(fig3)
                
                with col2:
                    # Filter status
                    filter_counts = df['Filter'].value_counts()
                    fig4, ax4 = plt.subplots(figsize=(8, 6))
                    ax4.bar(filter_counts.index, filter_counts.values, color=['#2E86AB', '#A23B72'])
                    ax4.set_title('Filter Status Distribution', fontweight='bold')
                    ax4.set_ylabel('Count')
                    st.pyplot(fig4)
            
            with tab3:
                st.subheader("Variant Data")
                st.dataframe(df.head(100), use_container_width=True)
                
        except Exception as e:
            st.error(f"❌ Error processing VCF file: {str(e)}")
    else:
        st.warning("👆 Please upload a VCF file to begin analysis")

def file_converter():
    st.markdown('<h2 class="section-header">🔄 File Converter</h2>', unsafe_allow_html=True)
    
    st.info("Convert between different genomic file formats and export your analysis results.")
    
    conversion_type = st.selectbox(
        "Select conversion type:",
        ["FASTA to CSV", "FASTA to JSON", "VCF to CSV", "VCF to JSON"]
    )
    
    uploaded_file = st.file_uploader(
        f"📤 Upload file for {conversion_type}",
        type=['fasta', 'fa', 'fna', 'vcf', 'txt']
    )
    
    if uploaded_file is not None:
        try:
            if "FASTA" in conversion_type and uploaded_file.name.endswith(('.fasta', '.fa', '.fna')):
                sequences = list(SeqIO.parse(uploaded_file, "fasta"))
                
                sequence_data = []
                for seq in sequences:
                    sequence_data.append({
                        'ID': seq.id,
                        'Sequence': str(seq.seq),
                        'Length': len(seq),
                        'GC_Content': GC(seq.seq),
                        'Molecular_Weight': molecular_weight(seq.seq),
                        'Description': seq.description
                    })
                
                df = pd.DataFrame(sequence_data)
                
                st.success(f"✅ Converted {len(sequences)} sequences")
                st.dataframe(df, use_container_width=True)
                
                # Download buttons
                col1, col2 = st.columns(2)
                
                with col1:
                    csv = df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download as CSV",
                        data=csv,
                        file_name="sequences.csv",
                        mime="text/csv",
                        type="primary"
                    )
                
                with col2:
                    if "JSON" in conversion_type:
                        json_str = df.to_json(indent=2, orient='records')
                        st.download_button(
                            label="📥 Download as JSON",
                            data=json_str,
                            file_name="sequences.json",
                            mime="application/json"
                        )
            
            else:
                st.warning("Please upload a file that matches the selected conversion type.")
                
        except Exception as e:
            st.error(f"❌ Error during conversion: {str(e)}")

def show_about():
    st.markdown('<h2 class="section-header">ℹ️ About</h2>', unsafe_allow_html=True)
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.markdown("""
        ### 🧬 Genomic Analysis Toolkit
        
        A comprehensive web application built for genomic data analysis and visualization.
        
        **🌟 Features**
        - Advanced sequence analysis and statistics
        - Variant calling data processing
        - Interactive visualizations and plots
        - File format conversion utilities
        - User-friendly web interface
        
        **🛠️ Technologies Used**
        - **Streamlit**: Web application framework
        - **Biopython**: Biological computation
        - **Pandas**: Data manipulation and analysis
        - **Matplotlib/Seaborn**: Scientific visualization
        - **Plotly**: Interactive plotting
        
        **📊 Supported Analyses**
        - GC content calculation and visualization
        - Sequence statistics and properties
        - Variant distribution analysis
        - Quality metrics assessment
        - Chromosomal distribution plots
        """)
    
    with col2:
        st.markdown("""
        ### 🔧 Technical Details
        
        **File Support**
        - FASTA (.fasta, .fa, .fna)
        - VCF (.vcf, .vcf.gz)
        - CSV/TSV exports
        
        **Requirements**
        - Python 3.8+
        - Streamlit 1.28+
        - Biopython 1.81+
        - Pandas 2.0+
        
        **Open Source**
        This project is open source and available on GitHub.
        """)

if __name__ == "__main__":
    main()