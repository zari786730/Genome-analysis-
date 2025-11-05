from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Dict, Optional
import uvicorn
import io
import json
from datetime import datetime

# Initialize FastAPI app
app = FastAPI(
    title="Genomic Analysis API",
    description="A comprehensive DNA/RNA sequence analysis API",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Request models
class SequenceAnalysisRequest(BaseModel):
    sequence: str
    analyses: List[str]
    sequence_type: str = "DNA"
    circular: bool = False

class FileConversionRequest(BaseModel):
    input_format: str
    output_format: str
    sequences: List[str]

class AdvancedAnalysisRequest(BaseModel):
    sequence: str
    analysis_type: str
    parameters: Optional[Dict] = {}

# Custom genomic analysis functions (no external dependencies)
def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content percentage"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total_count = len(sequence)
    
    if total_count == 0:
        return 0.0
    return round((gc_count / total_count) * 100, 2)

def calculate_molecular_weight(sequence: str, seq_type: str = "DNA") -> float:
    """Calculate molecular weight for DNA or RNA sequences"""
    sequence = sequence.upper()
    
    if seq_type == "DNA":
        # Molecular weights for DNA nucleotides (g/mol)
        weights = {'A': 331.2, 'T': 322.2, 'C': 307.2, 'G': 347.2}
    else:  # RNA
        weights = {'A': 347.2, 'U': 324.2, 'C': 323.2, 'G': 363.2}
    
    total_weight = 0.0
    valid_nucleotides = 0
    
    for nucleotide in sequence:
        if nucleotide in weights:
            total_weight += weights[nucleotide]
            valid_nucleotides += 1
    
    # Add water molecular weight for the chain (18 g/mol per bond)
    if valid_nucleotides > 1:
        total_weight += 18.0 * (valid_nucleotides - 1)
    
    return round(total_weight, 2)

def get_reverse_complement(sequence: str, seq_type: str = "DNA") -> str:
    """Get reverse complement of DNA/RNA sequence"""
    sequence = sequence.upper()
    
    if seq_type == "DNA":
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    else:  # RNA
        complement_map = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
    
    complement = ''
    for nucleotide in sequence:
        complement += complement_map.get(nucleotide, nucleotide)
    
    return complement[::-1]  # Reverse the complement

def calculate_melting_temperature(sequence: str) -> float:
    """Calculate melting temperature using Wallace rule"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return round(2 * gc_count + 4 * (len(sequence) - gc_count), 2)

def translate_dna_to_protein(sequence: str) -> str:
    """Translate DNA sequence to protein"""
    sequence = sequence.upper()
    
    # Standard genetic code
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
    }
    
    protein = ''
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            protein += codon_table.get(codon, 'X')
    
    return protein

def find_restriction_sites(sequence: str) -> Dict:
    """Find common restriction enzyme cutting sites"""
    sequence = sequence.upper()
    
    restriction_enzymes = {
        'EcoRI': 'GAATTC',
        'BamHI': 'GGATCC',
        'HindIII': 'AAGCTT',
        'XbaI': 'TCTAGA',
        'NotI': 'GCGGCCGC',
        'SacI': 'GAGCTC',
        'KpnI': 'GGTACC',
        'PstI': 'CTGCAG',
        'SmaI': 'CCCGGG'
    }
    
    sites_found = {}
    
    for enzyme, site in restriction_enzymes.items():
        positions = []
        start = 0
        while start < len(sequence):
            pos = sequence.find(site, start)
            if pos == -1:
                break
            positions.append(pos + 1)  # Convert to 1-based indexing
            start = pos + 1
        
        if positions:
            sites_found[enzyme] = positions
    
    return sites_found

def find_orfs(sequence: str, min_length: int = 50) -> List[Dict]:
    """Find Open Reading Frames"""
    sequence = sequence.upper()
    orfs = []
    
    # Check all reading frames
    for frame in range(3):
        # Forward strand
        protein = translate_dna_to_protein(sequence[frame:])
        orfs.extend(_find_orfs_in_translation(protein, frame + 1, "forward", min_length))
        
        # Reverse complement strand
        rev_comp = get_reverse_complement(sequence, "DNA")
        rev_protein = translate_dna_to_protein(rev_comp[frame:])
        orfs.extend(_find_orfs_in_translation(rev_protein, frame + 1, "reverse", min_length))
    
    # Sort by length and return top 10
    orfs.sort(key=lambda x: x['length'], reverse=True)
    return orfs[:10]

def _find_orfs_in_translation(protein: str, frame: int, strand: str, min_length: int) -> List[Dict]:
    """Helper function to find ORFs in translated sequence"""
    orfs = []
    start_pos = 0
    
    while start_pos < len(protein):
        # Find start codon (M)
        start_index = protein.find('M', start_pos)
        if start_index == -1:
            break
        
        # Find stop codon (*) after start
        stop_index = protein.find('*', start_index)
        if stop_index == -1:
            stop_index = len(protein)
        
        orf_length = (stop_index - start_index) * 3
        if orf_length >= min_length:
            orfs.append({
                "start": frame + (start_index * 3),
                "length": orf_length,
                "sequence": protein[start_index:stop_index],
                "strand": strand
            })
        
        start_pos = stop_index + 1
    
    return orfs

def calculate_gc_skew(sequence: str, window: int = 100) -> List[float]:
    """Calculate GC skew across the sequence"""
    sequence = sequence.upper()
    gc_skew = []
    
    for i in range(0, len(sequence) - window + 1, window):
        window_seq = sequence[i:i + window]
        g_count = window_seq.count('G')
        c_count = window_seq.count('C')
        
        if g_count + c_count > 0:
            skew = (g_count - c_count) / (g_count + c_count)
        else:
            skew = 0.0
        
        gc_skew.append(round(skew, 3))
    
    return gc_skew

def calculate_codon_usage(sequence: str) -> Dict[str, float]:
    """Calculate codon usage frequencies"""
    sequence = sequence.upper()
    codon_counts = {}
    total_codons = 0
    
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if len(codon) == 3 and all(nuc in 'ATCG' for nuc in codon):
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
            total_codons += 1
    
    if total_codons > 0:
        codon_freq = {codon: round(count/total_codons, 4) for codon, count in codon_counts.items()}
    else:
        codon_freq = {}
    
    return codon_freq

def parse_fasta_content(content: str) -> List[Dict]:
    """Parse FASTA format content"""
    sequences = []
    current_seq = ""
    current_id = ""
    current_desc = ""
    
    lines = content.split('\n')
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence
            if current_seq and current_id:
                sequences.append({
                    "id": current_id,
                    "description": current_desc,
                    "sequence": current_seq,
                    "length": len(current_seq)
                })
            
            # Start new sequence
            header = line[1:].strip()
            parts = header.split(' ', 1)
            current_id = parts[0]
            current_desc = parts[1] if len(parts) > 1 else ""
            current_seq = ""
        else:
            current_seq += line
    
    # Add the last sequence
    if current_seq and current_id:
        sequences.append({
            "id": current_id,
            "description": current_desc,
            "sequence": current_seq,
            "length": len(current_seq)
        })
    
    return sequences

# API Routes
@app.get("/")
async def root():
    return {
        "message": "Genomic Analysis API",
        "version": "1.0.0",
        "status": "running",
        "timestamp": datetime.now().isoformat()
    }

@app.get("/health")
async def health_check():
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}

@app.post("/analyze/sequence")
async def analyze_sequence(request: SequenceAnalysisRequest):
    """Analyze DNA/RNA sequences for various properties"""
    try:
        sequence = request.sequence.upper().strip()
        if not all(c in 'ATCGNU ' for c in sequence):
            raise HTTPException(status_code=400, detail="Invalid sequence characters")
        
        clean_sequence = sequence.replace(' ', '')
        if len(clean_sequence) < 10:
            raise HTTPException(status_code=400, detail="Sequence too short (minimum 10 bp)")
        
        results = {}
        
        for analysis in request.analyses:
            if analysis == "GC Content":
                results['gc_content'] = calculate_gc_content(clean_sequence)
            
            elif analysis == "Molecular Weight":
                results['molecular_weight'] = calculate_molecular_weight(clean_sequence, request.sequence_type)
            
            elif analysis == "Sequence Length":
                results['sequence_length'] = len(clean_sequence)
            
            elif analysis == "Nucleotide Frequency":
                results['nucleotide_frequency'] = {
                    'A': clean_sequence.count('A'),
                    'T': clean_sequence.count('T'),
                    'C': clean_sequence.count('C'),
                    'G': clean_sequence.count('G'),
                    'U': clean_sequence.count('U') if request.sequence_type == "RNA" else 0
                }
            
            elif analysis == "Reverse Complement":
                results['reverse_complement'] = get_reverse_complement(clean_sequence, request.sequence_type)
            
            elif analysis == "Melting Temperature":
                results['melting_temperature'] = calculate_melting_temperature(clean_sequence)
            
            elif analysis == "Translation":
                if request.sequence_type == "DNA":
                    results['translation'] = translate_dna_to_protein(clean_sequence)
            
            elif analysis == "Restriction Sites":
                results['restriction_sites'] = find_restriction_sites(clean_sequence)
        
        return JSONResponse(content=results)
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis error: {str(e)}")

@app.post("/convert/format")
async def convert_sequence_format(request: FileConversionRequest):
    """Convert between sequence formats"""
    try:
        converted_sequences = []
        
        for seq in request.sequences:
            if request.input_format.upper() == "FASTA" and request.output_format.upper() == "GENBANK":
                try:
                    # Parse FASTA
                    sequences = parse_fasta_content(seq)
                    if sequences:
                        record = sequences[0]
                        converted = f"LOCUS       {record['id']} {record['length']} bp DNA\n"
                        converted += f"DEFINITION  {record['description']}\n"
                        converted += f"ORIGIN\n"
                        
                        seq_str = record['sequence']
                        for i in range(0, len(seq_str), 60):
                            block = seq_str[i:i+60]
                            line_num = i + 1
                            converted += f"{line_num:9} {block}\n"
                        
                        converted += "//"
                        converted_sequences.append(converted)
                except Exception as e:
                    converted_sequences.append(f"Error converting sequence: {str(e)}")
            
            elif request.input_format.upper() == "GENBANK" and request.output_format.upper() == "FASTA":
                # Simple conversion - just extract sequence and create FASTA
                try:
                    lines = seq.split('\n')
                    sequence_data = ""
                    in_origin = False
                    
                    for line in lines:
                        if line.startswith("ORIGIN"):
                            in_origin = True
                            continue
                        elif in_origin and line.startswith("//"):
                            break
                        elif in_origin:
                            # Remove numbers and spaces
                            sequence_part = ''.join([c for c in line if c in 'atcgATCG'])
                            sequence_data += sequence_part
                    
                    converted = f">converted_sequence\n{sequence_data}"
                    converted_sequences.append(converted)
                except Exception as e:
                    converted_sequences.append(f"Error converting sequence: {str(e)}")
            else:
                converted_sequences.append(f"Unsupported conversion: {request.input_format} to {request.output_format}")
        
        return {"converted_sequences": converted_sequences}
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Conversion error: {str(e)}")

@app.post("/analyze/advanced")
async def advanced_analysis(request: AdvancedAnalysisRequest):
    """Perform advanced genomic analyses"""
    try:
        sequence = request.sequence.upper().strip()
        clean_sequence = sequence.replace(' ', '')
        
        if len(clean_sequence) < 10:
            raise HTTPException(status_code=400, detail="Sequence too short for advanced analysis")
        
        results = {}
        
        if request.analysis_type == "orf_finder":
            min_length = request.parameters.get("min_length", 50) if request.parameters else 50
            results['orfs'] = find_orfs(clean_sequence, min_length)
        
        elif request.analysis_type == "gc_skew":
            results['gc_skew'] = calculate_gc_skew(clean_sequence)
        
        elif request.analysis_type == "codon_usage":
            results['codon_usage'] = calculate_codon_usage(clean_sequence)
        else:
            raise HTTPException(status_code=400, detail="Unsupported analysis type")
        
        return JSONResponse(content=results)
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Advanced analysis error: {str(e)}")

@app.post("/upload/file")
async def upload_sequence_file(file: UploadFile = File(...)):
    """Upload and process sequence files"""
    try:
        content = await file.read()
        file_type = file.filename.split('.')[-1].lower()
        
        sequences = []
        
        if file_type in ['fasta', 'fa', 'fna']:
            sequences = parse_fasta_content(content.decode())
        
        else:
            raise HTTPException(status_code=400, detail="Unsupported file format. Only FASTA files are supported.")
        
        return {
            "filename": file.filename,
            "file_type": file_type,
            "sequences": sequences,
            "total_sequences": len(sequences)
        }
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"File processing error: {str(e)}")

# Main execution
if __name__ == "__main__":
    uvicorn.run(
        "backend:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )