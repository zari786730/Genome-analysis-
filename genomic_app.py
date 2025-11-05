import streamlit as st
import requests
import json
import plotly.graph_objects as go
import plotly.express as px
from datetime import datetime
import base64

# Page configuration
st.set_page_config(
    page_title="Genomic Analysis Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for styling
def local_css():
    st.markdown("""
    <style>
    .main-header {
        font-size: 3rem;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 2rem;
    }
    .feature-card {
        background-color: #F8F9FA;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 4px solid #2E86AB;
        margin: 1rem 0;
    }
    .result-box {
        background-color: #E8F4F8;
        padding: 1rem;
        border-radius: 8px;
        border: 1px solid #2E86AB;
    }
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 1rem;
        border-radius: 10px;
        text-align: center;
    }
    .stButton>button {
        width: 100%;
    }
    </style>
    """, unsafe_allow_html=True)

local_css()

# Backend API URL
API_URL = "http://localhost:8000"

def main():
    # Header Section
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        st.markdown('<h1 class="main-header">🧬 Genomic Analysis Pro</h1>', unsafe_allow_html=True)
        st.markdown("### Professional DNA/RNA Sequence Analysis Platform")
    
    # Sidebar Navigation
    st.sidebar.image("https://img.icons8.com/color/96/000000/dna-helix.png", width=80)
    st.sidebar.title("Navigation")
    app_mode = st.sidebar.selectbox(
        "Choose Analysis Type",
        ["🏠 Dashboard", "🔍 Sequence Analysis", "📊 Advanced Analytics", 
         "🔄 File Converter", "📈 Visualization", "⚙️ API Documentation"]
    )
    
    if app_mode == "🏠 Dashboard":
        show_dashboard()
    elif app_mode == "🔍 Sequence Analysis":
        sequence_analysis()
    elif app_mode == "📊 Advanced Analytics":
        advanced_analytics()
    elif app_mode == "🔄 File Converter":
        file_converter()
    elif app_mode == "📈 Visualization":
        visualization()
    elif app_mode == "⚙️ API Documentation":
        api_documentation()

def show_dashboard():
    st.markdown("## 📊 Dashboard Overview")
    
    # Quick stats
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown('<div class="metric-card"><h3>GC Analysis</h3><h2>Ready</h2></div>', unsafe_allow_html=True)
    with col2:
        st.markdown('<div class="metric-card"><h3>MW Calculator</h3><h2>Ready</h2></div>', unsafe_allow_html=True)
    with col3:
        st.markdown('<div class="metric-card"><h3>File Converter</h3><h2>Ready</h2></div>', unsafe_allow_html=True)
    with col4:
        st.markdown('<div class="metric-card"><h3>Visualization</h3><h2>Ready</h2></div>', unsafe_allow_html=True)
    
    # Features section
    st.markdown("## 🚀 Features")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("""
        <div class="feature-card">
        <h4>🔍 Basic Sequence Analysis</h4>
        <ul>
        <li>GC Content Calculation</li>
        <li>Molecular Weight</li>
        <li>Sequence Length</li>
        <li>Nucleotide Frequency</li>
        <li>Reverse Complement</li>
        <li>Melting Temperature</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="feature-card">
        <h4>🔄 File Conversion</h4>
        <ul>
        <li>FASTA Format Support</li>
        <li>GenBank Conversion</li>
        <li>Multi-sequence Files</li>
        <li>Batch Processing</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="feature-card">
        <h4>📊 Advanced Analytics</h4>
        <ul>
        <li>Restriction Site Analysis</li>
        <li>ORF Finder</li>
        <li>Translation Tools</li>
        <li>Codon Usage Analysis</li>
        <li>GC Skew Analysis</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("""
        <div class="feature-card">
        <h4>📈 Visualization</h4>
        <ul>
        <li>Interactive Plots</li>
        <li>Sequence Composition</li>
        <li>GC Skew Analysis</li>
        <li>Custom Charts</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # Quick start example
    st.markdown("## ⚡ Quick Start")
    example_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG"
    if st.button("Try Example Sequence", use_container_width=True):
        st.session_state.quick_start = example_sequence
        st.rerun()

def sequence_analysis():
    st.markdown("## 🔍 Sequence Analysis")
    
    tab1, tab2 = st.tabs(["📝 Input Sequence", "📁 Upload File"])
    
    with tab1:
        col1, col2 = st.columns([2, 1])
        
        with col1:
            # Check if quick start sequence exists
            default_sequence = st.session_state.get('quick_start', '')
            
            sequence = st.text_area(
                "Enter DNA/RNA Sequence",
                height=200,
                value=default_sequence,
                placeholder="ATCGATCGATCGATCG...",
                help="Enter nucleotide sequence (A, T, C, G for DNA; A, U, C, G for RNA)"
            )
            
            analysis_options = st.multiselect(
                "Select Analyses",
                ["GC Content", "Molecular Weight", "Sequence Length", "Nucleotide Frequency", 
                 "Reverse Complement", "Melting Temperature", "Translation", "Restriction Sites"],
                default=["GC Content", "Molecular Weight", "Sequence Length", "Nucleotide Frequency"]
            )
        
        with col2:
            st.markdown("### 🎛️ Settings")
            sequence_type = st.radio("Sequence Type", ["DNA", "RNA"])
            circular = st.checkbox("Circular Sequence")
            
            st.markdown("### 💡 Tips")
            st.info("""
            - Use A, T, C, G for DNA sequences
            - Use A, U, C, G for RNA sequences  
            - Remove spaces and line breaks for best results
            - Minimum sequence length: 10 bp
            """)
            
            if st.button("🚀 Analyze Sequence", type="primary", use_container_width=True):
                if sequence and len(sequence.strip()) >= 10:
                    with st.spinner("Analyzing sequence..."):
                        try:
                            # Clean the sequence
                            clean_sequence = ''.join(sequence.upper().split())
                            
                            response = requests.post(
                                f"{API_URL}/analyze/sequence",
                                json={
                                    "sequence": clean_sequence,
                                    "analyses": analysis_options,
                                    "sequence_type": sequence_type,
                                    "circular": circular
                                },
                                timeout=30
                            )
                            
                            if response.status_code == 200:
                                results = response.json()
                                display_sequence_results(results)
                            else:
                                st.error(f"Analysis failed: {response.text}")
                        except requests.exceptions.ConnectionError:
                            st.error("❌ Cannot connect to backend server. Make sure the backend is running on port 8000.")
                        except Exception as e:
                            st.error(f"Analysis error: {str(e)}")
                else:
                    st.warning("Please enter a valid sequence (minimum 10 characters)")

    with tab2:
        st.markdown("### 📁 Upload Sequence File")
        uploaded_file = st.file_uploader("Choose a FASTA file", type=['fasta', 'fa', 'fna'])
        
        if uploaded_file is not None:
            try:
                # Read and display file content
                content = uploaded_file.getvalue().decode()
                st.text_area("File Content", content, height=150)
                
                if st.button("Analyze Uploaded File"):
                    with st.spinner("Processing file..."):
                        files = {"file": (uploaded_file.name, content, "text/plain")}
                        response = requests.post(f"{API_URL}/upload/file", files=files)
                        
                        if response.status_code == 200:
                            file_results = response.json()
                            st.success(f"✅ Successfully processed {len(file_results['sequences'])} sequences")
                            
                            for seq in file_results['sequences']:
                                with st.expander(f"Sequence: {seq['id']}"):
                                    st.write(f"Length: {seq['length']} bp")
                                    st.code(seq['sequence'][:100] + "..." if len(seq['sequence']) > 100 else seq['sequence'])
            except Exception as e:
                st.error(f"File processing error: {str(e)}")

def display_sequence_results(results):
    st.markdown("## 📊 Analysis Results")
    
    # Basic metrics in cards
    col1, col2, col3, col4 = st.columns(4)
    
    if 'gc_content' in results:
        with col1:
            st.metric("GC Content", f"{results['gc_content']:.2f}%")
    
    if 'molecular_weight' in results:
        with col2:
            st.metric("Molecular Weight", f"{results['molecular_weight']:.2f} g/mol")
    
    if 'sequence_length' in results:
        with col3:
            st.metric("Sequence Length", f"{results['sequence_length']} bp")
    
    if 'melting_temperature' in results:
        with col4:
            st.metric("Tm", f"{results['melting_temperature']:.2f}°C")
    
    # Detailed results
    if 'nucleotide_frequency' in results:
        st.markdown("### 🧬 Nucleotide Composition")
        freq = results['nucleotide_frequency']
        
        # Create two columns for chart and table
        col1, col2 = st.columns([2, 1])
        
        with col1:
            # Create pie chart
            labels = list(freq.keys())
            values = list(freq.values())
            
            fig = px.pie(
                values=values,
                names=labels,
                title="Nucleotide Distribution",
                color=labels,
                color_discrete_map={'A': '#FF6B6B', 'T': '#4ECDC4', 'C': '#45B7D1', 'G': '#96CEB4', 'U': '#FFEAA7'}
            )
            st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            # Display frequency table
            st.markdown("**Frequency Table:**")
            for nt, count in freq.items():
                if count > 0:
                    percentage = (count / sum(freq.values())) * 100
                    st.write(f"{nt}: {count} ({percentage:.1f}%)")
    
    if 'reverse_complement' in results:
        st.markdown("### 🔁 Reverse Complement")
        st.code(results['reverse_complement'])
    
    if 'translation' in results:
        st.markdown("### 🧪 Translation")
        st.code(results['translation'])
    
    if 'restriction_sites' in results:
        st.markdown("### ✂️ Restriction Sites")
        sites = results['restriction_sites']
        if sites:
            for enzyme, positions in sites.items():
                st.write(f"**{enzyme}**: {positions}")
        else:
            st.info("No restriction sites found for common enzymes")

def advanced_analytics():
    st.markdown("## 📊 Advanced Analytics")
    
    tab1, tab2, tab3 = st.tabs(["ORF Finder", "GC Skew Analysis", "Codon Usage"])
    
    with tab1:
        st.markdown("### Open Reading Frame (ORF) Finder")
        orf_sequence = st.text_area("Enter sequence for ORF analysis", height=150)
        min_length = st.slider("Minimum ORF length (bp)", 30, 300, 50)
        
        if st.button("Find ORFs"):
            if orf_sequence:
                try:
                    response = requests.post(
                        f"{API_URL}/analyze/advanced",
                        json={
                            "sequence": orf_sequence,
                            "analysis_type": "orf_finder",
                            "parameters": {"min_length": min_length}
                        }
                    )
                    if response.status_code == 200:
                        results = response.json()
                        display_orf_results(results)
                    else:
                        st.error("ORF analysis failed")
                except Exception as e:
                    st.error(f"ORF analysis error: {str(e)}")
    
    with tab2:
        st.markdown("### GC Skew Analysis")
        st.info("GC Skew analysis coming in next version!")
    
    with tab3:
        st.markdown("### Codon Usage Analysis")
        st.info("Codon usage analysis coming in next version!")

def display_orf_results(results):
    if 'orfs' in results:
        orfs = results['orfs']
        st.write(f"Found {len(orfs)} ORFs")
        
        for i, orf in enumerate(orfs, 1):
            with st.expander(f"ORF {i}: {orf['length']} bp (Strand: {orf['strand']})"):
                st.write(f"Start position: {orf['start']}")
                st.write(f"Length: {orf['length']} bp")
                st.write(f"Amino acid sequence: {orf['sequence']}")

def file_converter():
    st.markdown("## 🔄 File Converter")
    
    col1, col2 = st.columns(2)
    
    with col1:
        input_format = st.selectbox("Input Format", ["FASTA", "GenBank"])
        input_sequences = st.text_area("Input Sequences", height=200)
    
    with col2:
        output_format = st.selectbox("Output Format", ["GenBank", "FASTA"])
        
        if st.button("Convert Format"):
            if input_sequences:
                try:
                    sequences = input_sequences.split('>')[1:]  # Split FASTA sequences
                    sequences = [f">{seq}" for seq in sequences]
                    
                    response = requests.post(
                        f"{API_URL}/convert/format",
                        json={
                            "input_format": input_format,
                            "output_format": output_format,
                            "sequences": sequences
                        }
                    )
                    
                    if response.status_code == 200:
                        results = response.json()
                        st.success("Conversion successful!")
                        
                        for i, converted in enumerate(results['converted_sequences']):
                            st.text_area(f"Converted Sequence {i+1}", converted, height=150)
                    else:
                        st.error("Conversion failed")
                except Exception as e:
                    st.error(f"Conversion error: {str(e)}")

def visualization():
    st.markdown("## 📈 Visualization")
    st.info("Advanced visualization features coming in the next version!")
    
    # Placeholder for future visualizations
    st.write("Future features will include:")
    st.write("• Interactive sequence maps")
    st.write("• GC content sliding window plots")
    st.write("• Restriction site mapping")
    st.write("• Phylogenetic tree visualization")

def api_documentation():
    st.markdown("## ⚙️ API Documentation")
    
    st.markdown("""
    ### Backend API Endpoints
    
    Our genomic analysis platform provides a RESTful API with the following endpoints:
    
    #### 🔍 Sequence Analysis
    - `POST /analyze/sequence` - Analyze DNA/RNA sequences
    - `POST /analyze/advanced` - Advanced analyses (ORF finder, GC skew, etc.)
    
    #### 🔄 File Operations
    - `POST /upload/file` - Upload and process sequence files
    - `POST /convert/format` - Convert between file formats
    
    #### ℹ️ System Info
    - `GET /` - API information
    - `GET /health` - Health check
    
    ### Interactive API Docs
    Visit the interactive API documentation at: 
    **http://localhost:8000/docs**
    
    ### Example Usage
    ```python
    import requests
    
    response = requests.post(
        "http://localhost:8000/analyze/sequence",
        json={
            "sequence": "ATCGATCGATCG",
            "analyses": ["GC Content", "Molecular Weight"],
            "sequence_type": "DNA"
        }
    )
    results = response.json()
    ```
    """)

if __name__ == "__main__":
    main()