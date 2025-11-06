from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Dict, Optional
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

# Custom genomic analysis functions
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
        weights = {'A': 331.2, 'T': 322.2, 'C': 307.2, 'G': 347.2}
    else:  # RNA
        weights = {'A': 347.2, 'U': 324.2, 'C': 323.2, 'G': 363.2}
    
    total_weight = 0.0
    valid_nucleotides = 0
    
    for nucleotide in sequence:
        if nucleotide in weights:
            total_weight += weights[nucleotide]
            valid_nucleotides += 1
    
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
    
    return complement[::-1]

def calculate_melting_temperature(sequence: str) -> float:
    """Calculate melting temperature using Wallace rule"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return round(2 * gc_count + 4 * (len(sequence) - gc_count), 2)

def translate_dna_to_protein(sequence: str) -> str:
    """Translate DNA sequence to protein"""
    sequence = sequence.upper()
    
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
        'XbaI': 'TCTAGA'
    }
    
    sites_found = {}
    
    for enzyme, site in restriction_enzymes.items():
        positions = []
        start = 0
        while start < len(sequence):
            pos = sequence.find(site, start)
            if pos == -1:
                break
            positions.append(pos + 1)
            start = pos + 1
        
        if positions:
            sites_found[enzyme] = positions
    
    return sites_found

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
            if current_seq and current_id:
                sequences.append({
                    "id": current_id,
                    "description": current_desc,
                    "sequence": current_seq,
                    "length": len(current_seq)
                })
            
            header = line[1:].strip()
            parts = header.split(' ', 1)
            current_id = parts[0]
            current_desc = parts[1] if len(parts) > 1 else ""
            current_seq = ""
        else:
            current_seq += line
    
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