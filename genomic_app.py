import streamlit as st
import requests
import plotly.express as px

# Page configuration
st.set_page_config(
    page_title="Genomic Analysis Pro",
    page_icon="🧬",
    layout="wide"
)

# Backend API URL
API_URL = "http://localhost:8000"

def main():
    st.title("🧬 Genomic Analysis Pro")
    st.markdown("### Professional DNA/RNA Sequence Analysis Platform")
    
    # Sidebar Navigation
    st.sidebar.title("Navigation")
    app_mode = st.sidebar.selectbox(
        "Choose Analysis Type",
        ["🏠 Dashboard", "🔍 Sequence Analysis", "📊 Advanced Analytics", "🔄 File Converter"]
    )
    
    if app_mode == "🏠 Dashboard":
        show_dashboard()
    elif app_mode == "🔍 Sequence Analysis":
        sequence_analysis()
    elif app_mode == "📊 Advanced Analytics":
        advanced_analytics()
    elif app_mode == "🔄 File Converter":
        file_converter()

def show_dashboard():
    st.markdown("## 📊 Dashboard Overview")
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("GC Analysis", "Ready")
    with col2:
        st.metric("MW Calculator", "Ready")
    with col3:
        st.metric("File Converter", "Ready")
    with col4:
        st.metric("ORF Finder", "Ready")
    
    st.markdown("## 🚀 Features")
    st.write("""
    - **GC Content Analysis**: Calculate GC percentage of sequences
    - **Molecular Weight**: Compute molecular weight of DNA/RNA
    - **Sequence Conversion**: Convert between FASTA and GenBank formats
    - **ORF Finding**: Identify open reading frames
    - **Restriction Sites**: Find enzyme cutting sites
    - **Codon Usage**: Analyze codon frequency patterns
    """)

def sequence_analysis():
    st.markdown("## 🔍 Sequence Analysis")
    
    sequence = st.text_area(
        "Enter DNA/RNA Sequence",
        height=200,
        placeholder="ATCGATCGATCGATCG...",
        help="Enter nucleotide sequence (A, T, C, G for DNA; A, U, C, G for RNA)"
    )
    
    col1, col2 = st.columns(2)
    
    with col1:
        analysis_options = st.multiselect(
            "Select Analyses",
            ["GC Content", "Molecular Weight", "Sequence Length", "Nucleotide Frequency", 
             "Reverse Complement", "Translation", "Restriction Sites"],
            default=["GC Content", "Molecular Weight", "Sequence Length"]
        )
    
    with col2:
        sequence_type = st.radio("Sequence Type", ["DNA", "RNA"])
        
        if st.button("🚀 Analyze Sequence", type="primary"):
            if sequence:
                with st.spinner("Analyzing sequence..."):
                    try:
                        clean_sequence = ''.join(sequence.upper().split())
                        
                        response = requests.post(
                            f"{API_URL}/analyze/sequence",
                            json={
                                "sequence": clean_sequence,
                                "analyses": analysis_options,
                                "sequence_type": sequence_type,
                                "circular": False
                            },
                            timeout=30
                        )
                        
                        if response.status_code == 200:
                            results = response.json()
                            display_sequence_results(results)
                        else:
                            st.error(f"Analysis failed: {response.text}")
                    except requests.exceptions.ConnectionError:
                        st.error("❌ Cannot connect to backend server. Make sure backend.py is running on port 8000.")
                    except Exception as e:
                        st.error(f"Analysis error: {str(e)}")
            else:
                st.warning("Please enter a sequence to analyze")

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
    
    # Nucleotide frequency chart
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
    
    tab1, tab2 = st.tabs(["ORF Finder", "File Upload"])
    
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
                        if 'orfs' in results:
                            orfs = results['orfs']
                            st.write(f"Found {len(orfs)} ORFs")
                            for i, orf in enumerate(orfs, 1):
                                with st.expander(f"ORF {i}: {orf['length']} bp"):
                                    st.write(f"Start: {orf['start']}, Strand: {orf['strand']}")
                                    st.code(orf['sequence'])
                    else:
                        st.error("ORF analysis failed")
                except Exception as e:
                    st.error(f"ORF analysis error: {str(e)}")
    
    with tab2:
        st.markdown("### 📁 Upload Sequence File")
        uploaded_file = st.file_uploader("Choose a FASTA file", type=['fasta', 'fa', 'fna'])
        
        if uploaded_file is not None:
            try:
                content = uploaded_file.getvalue().decode()
                st.text_area("File Content", content, height=150)
                
                if st.button("Process Uploaded File"):
                    with st.spinner("Processing file..."):
                        files = {"file": (uploaded_file.name, content, "text/plain")}
                        response = requests.post(f"{API_URL}/upload/file", files=files)
                        
                        if response.status_code == 200:
                            file_results = response.json()
                            st.success(f"✅ Successfully processed {file_results['total_sequences']} sequences")
            except Exception as e:
                st.error(f"File processing error: {str(e)}")

def file_converter():
    st.markdown("## 🔄 File Converter")
    
    col1, col2 = st.columns(2)
    
    with col1:
        input_format = st.selectbox("Input Format", ["FASTA", "GenBank"])
        input_sequences = st.text_area("Input Sequences", height=200, 
                                     placeholder=">sequence1\nATCG...\n>sequence2\nATCG...")
    
    with col2:
        output_format = st.selectbox("Output Format", ["GenBank", "FASTA"])
        
        if st.button("Convert Format"):
            if input_sequences:
                try:
                    sequences = [seq.strip() for seq in input_sequences.split('>') if seq.strip()]
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

if __name__ == "__main__":
    main()