
sequence = 'GTCCTGGCTGCTCTTCTAGAAGGTAATTCATGGAGAACAAGAGCTACTGAGGATTTGGCTGGTCGTGAGTGAGGGAATCAGAGGATGTGTGACAGTCTCCTGACCAGGATGTCTTTGTGTTTGCAGGTGTCCAGGCTCAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTCGGTGCAGGCTGGAGGGTCTCTGAAACTCTCCTGTGCAGCCTCTGGATACATCTTCAGTAGCTGCGGAATGGGCTGGTACCGCCAGGCTCCAGGGAAGGAGCGCGAGTTGGTCTCAACTATTAGTAGTGATGGTACCACAAGCTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCCAAGACAATGCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAAAACTGAGGACACGGCCATGTACTACTGT'



def translate(sequence):
    '''
    protein coding table:
    '''
    start_code = 'ATG'
    end_code = ['TAA', 'TAG', 'TGA']
    protein_table = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                     'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                     'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                     'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                     'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                     'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                     'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                     'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                     'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                     'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                     'TAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                     'TAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                     'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                     'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                     'TGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                     'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

    # Truncate sequence length to 3 * n
    seq = sequence[:len(sequence) // 3 * 3]

    protein = ''
    for i in range(0, len(seq), 3):
        protein = protein + protein_table[seq[i:i+3]]

    return protein


for i in range(0, 3):
    coding_sequence = sequence[130+i:]
    aa = translate(coding_sequence)
    # aa = str(Seq(coding_sequence, generic_dna).translate())

    # Has stop code in sequence
    has_stop = '*' in aa
    if not has_stop: print(aa)


