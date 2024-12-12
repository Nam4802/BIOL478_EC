from Bio import SeqIO
from Bio.Align import Alignment
from Bio.Seq import Seq
from algo import SW

def read_fasta(file_path):
    return  SeqIO.read(file_path, "fasta")

def dna2prot(dna_seq):
    return dna_seq.translate()

if __name__ == "__main__":
    dna_file = "dna_sequence.fasta"  # Replace with your input file
    protein_file = "protein_sequence.fasta"  # Replace with your input file
    gap_open_penalty = -7
    gap_extend_penalty = -1

    dna_seq = read_fasta(dna_file)
    prot_seq = read_fasta(protein_file)

    translated_prot_seq = dna2prot(dna_seq)

    # Running alignment with Smith - Waterman
    lines = SW(translated_prot_seq.seq, prot_seq.seq, gap_open_penalty, gap_extend_penalty)
    
    # Displaying alignment with Bio.Align.Alignment class
    sequences, coordinates = Alignment.parse_printed_alignment([bytes(''.join(line), 'utf-8') for line in lines])
    sequences = [Seq(sequence) for sequence in sequences]
    alignment = Alignment(sequences, coordinates)
    print(alignment)

