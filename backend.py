# app/backend.py
import pandas as pd
import numpy as np
from Bio import SeqIO, SeqUtils
from Bio.SeqUtils import GC, molecular_weight, MeltingTemp
from Bio.Seq import Seq
import cyvcf2
import tempfile
import os

class SequenceAnalyzer:
    def __init__(self):
        self.sequences = {}
    
    def load_sequences(self, file_content, file_type="fasta"):
        """Load sequences from uploaded file"""
        try:
            if file_type == "fasta":
                sequences = list(SeqIO.parse(file_content, "fasta"))
                for seq in sequences:
                    self.sequences[seq.id] = seq
                return True, f"Successfully loaded {len(sequences)} sequences"
            return False, "Unsupported file type"
        except Exception as e:
            return False, f"Error loading sequences: {str(e)}"
    
    def get_sequence_stats(self, sequence):
        """Get comprehensive statistics for a sequence"""
        stats = {
            'ID': sequence.id,
            'Length': len(sequence),
            'GC_Content': GC(sequence.seq),
            'Molecular_Weight': molecular_weight(sequence.seq),
            'Melting_Temp': MeltingTemp.Tm_Wallace(sequence.seq),
            'AT_Content': 100 - GC(sequence.seq)
        }
        
        # Nucleotide counts
        nucleotides = ['A', 'T', 'G', 'C']
        for nt in nucleotides:
            stats[f'count_{nt}'] = sequence.seq.count(nt)
        
        return stats
    
    def calculate_gc_content_sliding(self, sequence, window_size=100):
        """Calculate GC content in sliding windows"""
        gc_values = []
        positions = []
        
        seq_str = str(sequence.seq)
        for i in range(0, len(seq_str) - window_size + 1, window_size):
            window = seq_str[i:i + window_size]
            gc_percent = GC(window)
            gc_values.append(gc_percent)
            positions.append(i)
        
        return positions, gc_values
    
    def find_motifs(self, sequence, motif):
        """Find specific motifs in sequence"""
        positions = []
        motif = motif.upper()
        seq_str = str(sequence.seq).upper()
        
        start = 0
        while True:
            pos = seq_str.find(motif, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        
        return positions

class VariantAnalyzer:
    def __init__(self):
        self.variants_df = None
    
    def load_vcf_file(self, file_path):
        """Load and parse VCF file"""
        try:
            variants = []
            vcf = cyvcf2.VCF(file_path)
            
            for variant in vcf:
                variant_info = {
                    'CHROM': variant.CHROM,
                    'POS': variant.POS,
                    'ID': variant.ID if variant.ID else '.',
                    'REF': variant.REF,
                    'ALT': ','.join(variant.ALT),
                    'QUAL': variant.QUAL,
                    'FILTER': variant.FILTER,
                    'is_SNP': variant.is_snp,
                    'is_INDEL': variant.is_indel,
                    'is_Transition': variant.is_transition,
                    'num_Het': variant.num_het,
                    'num_Hom_Alt': variant.num_hom_alt,
                    'num_Hom_Ref': variant.num_hom_ref,
                    'call_Rate': variant.call_rate,
                    'AAF': variant.aaf
                }
                variants.append(variant_info)
            
            self.variants_df = pd.DataFrame(variants)
            return True, f"Loaded {len(self.variants_df)} variants"
            
        except Exception as e:
            return False, f"Error loading VCF file: {str(e)}"
    
    def get_variant_summary(self):
        """Get summary statistics for variants"""
        if self.variants_df is None:
            return None
        
        summary = {
            'Total Variants': len(self.variants_df),
            'SNPs': self.variants_df['is_SNP'].sum(),
            'INDELs': self.variants_df['is_INDEL'].sum(),
            'Transitions': self.variants_df['is_Transition'].sum(),
            'Average Quality': self.variants_df['QUAL'].mean(),
            'Average Call Rate': self.variants_df['call_Rate'].mean(),
            'Chromosomes': self.variants_df['CHROM'].nunique()
        }
        
        return summary
    
    def get_chromosome_distribution(self):
        """Get variant distribution by chromosome"""
        if self.variants_df is None:
            return None
        return self.variants_df['CHROM'].value_counts()

class GenomicUtilities:
    @staticmethod
    def complement_sequence(sequence):
        """Get complement of DNA sequence"""
        return str(Seq(sequence).complement())
    
    @staticmethod
    def reverse_complement(sequence):
        """Get reverse complement of DNA sequence"""
        return str(Seq(sequence).reverse_complement())
    
    @staticmethod
    def translate_sequence(sequence):
        """Translate DNA sequence to protein"""
        try:
            return str(Seq(sequence).translate())
        except:
            return "Translation failed - check sequence"

# Create singleton instances
sequence_analyzer = SequenceAnalyzer()
variant_analyzer = VariantAnalyzer()
genomic_utilities = GenomicUtilities()