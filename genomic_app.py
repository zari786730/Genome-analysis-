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
    </style>
    """, unsafe_allow_html=True)

local_css()

# Backend API URL (adjust for deployment)
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
        <li>Pattern Matching</li>
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

def sequence_analysis():
    st.markdown("## 🔍 Sequence Analysis")
    
    tab1, tab2, tab3 = st.tabs(["📝 Input Sequence", "📁 Upload File", "🎯 Batch Analysis"])
    
    with tab1:
        col1, col2 = st.columns([2, 1])
        
        with col1:
            sequence = st.text_area(
                "Enter DNA/RNA Sequence",
                height=200,
                placeholder="ATCGATCGATCGATCG...",
                help="Enter nucleotide sequence (A, T, C, G for DNA; A, U, C, G for RNA)"
            )
            
            analysis_options = st.multiselect(
                "Select Analyses",
                ["GC Content", "Molecular Weight", "Sequence Length", "Nucleotide Frequency", 
                 "Reverse Complement", "Translation", "Restriction Sites"],
                default=["GC Content", "Molecular Weight", "Sequence Length"]
            )
        
        with col2:
            st.markdown("### 🎛️ Settings")
            sequence_type = st.radio("Sequence Type", ["DNA", "RNA"])
            circular = st.checkbox("Circular Sequence")
            
            if st.button("🚀 Analyze Sequence", type="primary", use_container_width=True):
                if sequence:
                    with st.spinner("Analyzing sequence..."):
                        try:
                            response = requests.post(
                                f"{API_URL}/analyze/sequence",
                                json={
                                    "sequence": sequence,
                                    "analyses": analysis_options,
                                    "sequence_type": sequence_type,
                                    "circular": circular
                                }
                            )
                            
                            if response.status_code == 200:
                                results = response.json()
                                display_sequence_results(results)
                            else:
                                st.error("Analysis failed. Please check your sequence.")
                        except Exception as e:
                            st.error(f"Connection error: {e}")
                else:
                    st.warning("Please enter a sequence to analyze.")

def display_sequence_results(results):
    st.markdown("## 📊 Analysis Results")
    
    # Basic metrics
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
        fig = px.pie(
            values=list(freq.values()),
            names=list(freq.keys()),
            title="Nucleotide Distribution"
        )
        st.plotly_chart(fig, use_container_width=True)
    
    if 'reverse_complement' in results:
        st.markdown("### 🔁 Reverse Complement")
        st.code(results['reverse_complement'])

def advanced_analytics():
    st.markdown("## 📊 Advanced Analytics")
    # Implementation for advanced features
    st.info("Advanced analytics features coming soon!")

def file_converter():
    st.markdown("## 🔄 File Converter")
    # Implementation for file conversion
    st.info("File conversion features coming soon!")

def visualization():
    st.markdown("## 📈 Visualization")
    # Implementation for visualization
    st.info("Visualization features coming soon!")

def api_documentation():
    st.markdown("## ⚙️ API Documentation")
    st.info("API documentation available at /docs endpoint")

if __name__ == "__main__":
    main()