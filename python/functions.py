# Predefined functions for universal use


import pandas as pd




# Function to translate protein to DNA
def translate(sequence):
    '''
    protein coding table:
    '''
    start_code = 'ATG'
    end_code = ['TAA', 'TAG', 'TGA']
    protein_table = {'TTT': 'F', 
                     'CTT': 'L', 
                     'ATT': 'I', 
                     'GTT': 'V',
                     'TTC': 'F', 
                     'CTC': 'L', 
                     'ATC': 'I', 
                     'GTC': 'V',
                     'TTA': 'L', 
                     'CTA': 'L', 
                     'ATA': 'I', 
                     'GTA': 'V',
                     'TTG': 'L', 
                     'CTG': 'L', 
                     'ATG': 'M', 
                     'GTG': 'V',
                     'TCT': 'S', 
                     'CCT': 'P', 
                     'ACT': 'T', 
                     'GCT': 'A',
                     'TCC': 'S', 
                     'CCC': 'P', 
                     'ACC': 'T', 
                     'GCC': 'A',
                     'TCA': 'S', 
                     'CCA': 'P', 
                     'ACA': 'T', 
                     'GCA': 'A',
                     'TCG': 'S', 
                     'CCG': 'P', 
                     'ACG': 'T', 
                     'GCG': 'A',
                     'TAT': 'Y', 
                     'CAT': 'H', 
                     'AAT': 'N', 
                     'GAT': 'D',
                     'TAC': 'Y', 
                     'CAC': 'H', 
                     'AAC': 'N', 
                     'GAC': 'D',
                     'TAA': '*', 
                     'CAA': 'Q', 
                     'AAA': 'K', 
                     'GAA': 'E',
                     'TAG': '*', 
                     'CAG': 'Q', 
                     'AAG': 'K', 
                     'GAG': 'E',
                     'TGT': 'C', 
                     'CGT': 'R', 
                     'AGT': 'S', 
                     'GGT': 'G',
                     'TGC': 'C', 
                     'CGC': 'R', 
                     'AGC': 'S', 
                     'GGC': 'G',
                     'TGA': '*', 
                     'CGA': 'R', 
                     'AGA': 'R', 
                     'GGA': 'G',
                     'TGG': 'W', 
                     'CGG': 'R', 
                     'AGG': 'R', 
                     'GGG': 'G'}

    # Truncate sequence length to 3 * n
    seq = sequence[:len(sequence) // 3 * 3]

    protein = ''
    for i in range(0, len(seq), 3):
        protein = protein + protein_table[seq[i:i+3]]

    return protein


def nt2aa(nt_seq):
    '''
    Translate coding sequences of all 3 frames and select the longest as the final protein
    For the aa in each frame, extract the longest region without stop codon (between stop codon)
    '''
    # Tranlate coding sequnces in all 3 frames
    aa_seqs = [max(nt_seq[i:].translate().split("*"), key=len) for i in range(0,3)]
    
    # Select the longest frame as the final product
    aa_seq = max(aa_seqs, key=len)
   
    return(aa_seq)


def seq2matrix(seq, region_length):
    '''
    Split sequnce into subseqs (length = FR_length)
    Arguments: 
        seq: aa sequence
        region_length: FR length
    Returns: 
        Dataframe with subseqs in rows and aa in each column
    '''
    # Sequence length
    seq_length = len(seq)

    if seq_length == 0:
        return None

    # Convert string into letter list
    aa_list = [aa for aa in seq]

    # Initialize dataframe
    n_subseqs = abs(seq_length - region_length) + 1
    subseqs = pd.DataFrame(index=range(n_subseqs), columns=range(region_length), dtype=None)

    # Dump subseqs into dataframe
    if seq_length < region_length:
        for i in range(0, n_subseqs):
            subseq = [None] * i + aa_list + [None] * (region_length - seq_length - i)
            subseqs.loc[i] = subseq 
    elif seq_length >= region_length:
        for i in range(0, n_subseqs):
            subseq = aa_list[i: region_length + i]
            subseqs.loc[i] = subseq
    
    return subseqs


def get_heads(seq, region_len):
    '''
    Prepend gaps to the head
    Arguments:
      seq: query seq (str)
      region_len: FR region length (int)
    Return:
      subseqs (data.frame)
    '''
    n_inserts = int(0.3 * region_len)
    subseqs = pd.DataFrame(index=range(n_inserts), columns=range(region_len), dtype=None)
    for i in range(0, n_inserts):
        subseq = "-" * (i + 1) + seq[:region_len-i-1]
        subseq = [aa for aa in subseq]
        subseqs.loc[i] = subseq
    return(subseqs)


def get_tails(seq, region_len):
    ''''
    Append gaps to the tail
    Arguments:
      seq: query seq (str)
      region_len: RR region length (int)
    Return:
      subseqs (data.frame)
    '''
    n_inserts = int(0.3 * region_len)
    subseqs = pd.DataFrame(index=range(n_inserts), columns=range(region_len), dtype=None)
    for i in range(0, n_inserts):
        subseq = seq[-(region_len-i-1):] + "-" * (i + 1)
        subseq = [aa for aa in subseq]
        subseqs.loc[i] = subseq
    return(subseqs)






