from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Dict, Optional
import uvicorn
from Bio.SeqUtils import GC, molecular_weight, MeltingTemp
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Data import CodonTable
import io
import json
from datetime import datetime

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
    try:
        # Validate sequence
        sequence = request.sequence.upper().strip()
        if not all(c in 'ATCGNU ' for c in sequence):
            raise HTTPException(status_code=400, detail="Invalid sequence characters")
        
        seq_obj = Seq(sequence.replace(' ', ''))
        results = {}
        
        # Perform requested analyses
        for analysis in request.analyses:
            if analysis == "GC Content":
                results['gc_content'] = GC(seq_obj)
            
            elif analysis == "Molecular Weight":
                results['molecular_weight'] = molecular_weight(seq_obj, request.sequence_type == "DNA")
            
            elif analysis == "Sequence Length":
                results['sequence_length'] = len(seq_obj)
            
            elif analysis == "Nucleotide Frequency":
                results['nucleotide_frequency'] = {
                    'A': seq_obj.count('A'),
                    'T': seq_obj.count('T'),
                    'C': seq_obj.count('C'),
                    'G': seq_obj.count('G'),
                    'U': seq_obj.count('U') if request.sequence_type == "RNA" else 0
                }
            
            elif analysis == "Reverse Complement":
                results['reverse_complement'] = str(seq_obj.reverse_complement())
            
            elif analysis == "Melting Temperature":
                results['melting_temperature'] = MeltingTemp.Tm_Wallace(seq_obj)
            
            elif analysis == "Translation":
                if request.sequence_type == "DNA":
                    results['translation'] = str(seq_obj.translate())
            
            elif analysis == "Restriction Sites":
                results['restriction_sites'] = find_restriction_sites(seq_obj)
        
        return JSONResponse(content=results)
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis error: {str(e)}")

def find_restriction_sites(sequence: Seq) -> Dict:
    """Find restriction enzyme cutting sites"""
    try:
        rb = RestrictionBatch([], ['EcoRI', 'BamHI', 'HindIII', 'XbaI'])
        analysis = Analysis(rb, sequence)
        return {str(enzyme): sites for enzyme, sites in analysis.full().items() if sites}
    except:
        return {}

@app.post("/convert/format")
async def convert_sequence_format(request: FileConversionRequest):
    """Convert between sequence formats"""
    try:
        converted_sequences = []
        
        for seq in request.sequences:
            if request.input_format.upper() == "FASTA" and request.output_format.upper() == "GENBANK":
                # Convert FASTA to GenBank-like format
                record = SeqIO.read(io.StringIO(seq), "fasta")
                converted = f"LOCUS       {record.id} {len(record.seq)} bp DNA\n"
                converted += f"DEFINITION  {record.description}\n"
                converted += f"ORIGIN\n{record.seq}\n//"
                converted_sequences.append(converted)
            
            elif request.input_format.upper() == "GENBANK" and request.output_format.upper() == "FASTA":
                # Convert GenBank to FASTA
                record = SeqIO.read(io.StringIO(seq), "genbank")
                converted = f">{record.id} {record.description}\n{record.seq}"
                converted_sequences.append(converted)
        
        return {"converted_sequences": converted_sequences}
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Conversion error: {str(e)}")

@app.post("/analyze/advanced")
async def advanced_analysis(request: AdvancedAnalysisRequest):
    """Perform advanced genomic analyses"""
    try:
        sequence = request.sequence.upper().strip()
        seq_obj = Seq(sequence)
        results = {}
        
        if request.analysis_type == "orf_finder":
            results['orfs'] = find_orfs(seq_obj)
        
        elif request.analysis_type == "gc_skew":
            results['gc_skew'] = calculate_gc_skew(seq_obj)
        
        elif request.analysis_type == "codon_usage":
            results['codon_usage'] = calculate_codon_usage(seq_obj)
        
        return JSONResponse(content=results)
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Advanced analysis error: {str(e)}")

def find_orfs(sequence: Seq, min_length: int = 50) -> List[Dict]:
    """Find Open Reading Frames"""
    orfs = []
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(nuc) - frame) // 3)
            for pro in nuc[frame:frame+length].translate().split("*"):
                if len(pro) >= min_length // 3:
                    orfs.append({
                        "start": frame + 1,
                        "length": len(pro) * 3,
                        "sequence": str(pro),
                        "strand": "forward" if strand == 1 else "reverse"
                    })
    return orfs

def calculate_gc_skew(sequence: Seq, window: int = 100) -> List[float]:
    """Calculate GC skew across the sequence"""
    gc_skew = []
    for i in range(0, len(sequence) - window + 1, window):
        window_seq = sequence[i:i + window]
        g_count = window_seq.count('G')
        c_count = window_seq.count('C')
        if g_count + c_count > 0:
            skew = (g_count - c_count) / (g_count + c_count)
        else:
            skew = 0
        gc_skew.append(skew)
    return gc_skew

def calculate_codon_usage(sequence: Seq) -> Dict[str, float]:
    """Calculate codon usage frequencies"""
    codon_table = CodonTable.standard_dna_table
    codon_counts = {}
    total_codons = 0
    
    # Count codons
    for i in range(0, len(sequence) - 2, 3):
        codon = str(sequence[i:i+3])
        if len(codon) == 3:
            codon_counts[codon] = codon_counts.get(codon, 0) + 1
            total_codons += 1
    
    # Calculate frequencies
    codon_freq = {codon: count/total_codons for codon, count in codon_counts.items()}
    return codon_freq

@app.post("/upload/file")
async def upload_sequence_file(file: UploadFile = File(...)):
    """Upload and process sequence files"""
    try:
        content = await file.read()
        file_type = file.filename.split('.')[-1].lower()
        
        sequences = []
        if file_type in ['fasta', 'fa', 'fna']:
            records = SeqIO.parse(io.StringIO(content.decode()), "fasta")
            for record in records:
                sequences.append({
                    "id": record.id,
                    "description": record.description,
                    "sequence": str(record.seq),
                    "length": len(record.seq)
                })
        
        return {"filename": file.filename, "sequences": sequences}
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"File processing error: {str(e)}")

if __name__ == "__main__":
    uvicorn.run(
        "backend:app",
        host="0.0.0.0",
        port=8000,
        reload=True,
        log_level="info"
    )