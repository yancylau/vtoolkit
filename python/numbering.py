# IMGT unique numbering for FRs and CDRs
# Ref: http://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html#table1


def fr_numbering(seq, region):
    '''
    IMGT unique numbering for FRs
      FR1:   1 -  26
      FR2:  39 -  55
      FR3:  66 - 104
      FR4: 118 - 128
    Usually,
      * FR1 is 25 AAs with one gap in position 10;
      * FR2 is 12 AAs withnot gaps;
      * FR3 is 38 AAs with one gap in position 8; and
      * FR4 is 11 AAs without gaps
    Arguments: 
      seq: aa sequence (str) or None
      region: FR1/2/3/4 (str)
    Return:
      imgt numbered seq with gaps (str)
    '''
    # seq = seq.replace('-', '').replace('.', '')

    if region == "FR1":
        if len(seq) == 0: 
            return "." * 26
        if len(seq) == 25: 
            return seq[0:9] + "." + seq[9:25]
        if len(seq) == 26: 
            return seq
        return "X" * 26

    if region == "FR2":
        if len(seq) == 0:
            return "." * 17
        if len(seq) == 17:
            return seq
        return "X" * 17

    if region == "FR3":
        if len(seq) == 0:
            return "." * 39
        if len(seq) == 38:
            return seq[0:7] + "." + seq[7:38]
        if len(seq) == 39:
            return seq
        return "X" * 39

    if region == "FR4":
        if len(seq) == 0:
            return "." * 11
        if len(seq) == 11:
            return seq
        return "X" * 11   


def cdr_numbering(seq, region):
    '''
    IMGT unique numbering for CDRs
      CDR1          : 12  (5-12)
      CDR2          : 10  (0-10)
      CDR3 basic    : 13  (5-13)
      CDR3 extend 1 : 31 (14-31)
      CDR3 extend 2 : 51 (32-51)
      CDR3 extend 3 : 71 (52-71)
      CDR3 extend 4 : 91 (72-91)
    Arguments: 
      seq: aa sequence (str) or None
      region: CDR1/2 (str)
    Return:
      imgt numbered seq with gaps (str)
    '''
    # seq = seq.replace('-', '').replace('.', '')
    # seq_length = len(seq)

    # Gap orders in CDRs
    if region == 'CDR1':
        imgt_length = 12
        gaps = [7, 6, 8, 5, 9, 4, 10]

    if region == 'CDR2':
        imgt_length = 10
        gaps = [6, 5, 7, 4, 8, 3, 9, 2, 10, 1]

    if region == 'CDR3':
        imgt_length = 13
        gaps = [7, 8, 6, 9, 5, 10, 4, 11]
    if region == 'CDR3_extend_1':
        imgt_length = 31
        gaps = [16, 17, 15, 18, 14, 19, 13, 20, 12, 21, 11, 22, 10, 23, 9, 24, 8]
    if region == 'CDR3_extend_2':
        imgt_length = 51
        gaps = [26, 27, 25, 28, 24, 29, 23, 30, 22, 31, 21, 32, 20, 33, 19, 34, 18, 35, 17]
    if region == 'CDR3_extend_3':
        imgt_length = 71
        gaps = [36, 37, 35, 38, 34, 39, 33, 40, 32, 41, 31, 42, 30, 43, 29, 44, 28, 45, 27]
    if region == 'CDR3_extend_4':
        gaps = [46, 47, 45, 48, 44, 49, 43, 50, 42, 51, 41, 52, 40, 53, 39, 54, 38, 55, 37]
    
        
    # Add gaps to seq  
    if len(seq) == 0:
        imgt_seq = '.' * imgt_length
    elif len(seq) > imgt_length:
        # TODO:
        # May need to mark abnormal sequences that may arise from sequencing error.
        imgt_seq = seq
    else:
        # Find the positions to insert gaps
        n_gaps = imgt_length - int(len(seq))
        
        if n_gaps == 0: 
            imgt_seq =  seq
        else:
            gap_sites = gaps[:n_gaps]
            gaps_start = min(gap_sites) - 1
            # Insert gaps
            imgt_seq = seq[:gaps_start] + '.' * n_gaps + seq[gaps_start:]
    
    return imgt_seq



# IMGT unique numbering for FRs and CDRs
# Ref: http://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html#table1


def fr_numbering(seq, region):
    '''
    IMGT unique numbering for FRs
      FR1:   1 -  26
      FR2:  39 -  55
      FR3:  66 - 104
      FR4: 118 - 128
    Usually,
      * FR1 is 25 AAs with one gap in position 10;
      * FR2 is 12 AAs withnot gaps;
      * FR3 is 38 AAs with one gap in position 8; and
      * FR4 is 11 AAs without gaps
    Arguments: 
      seq: aa sequence (str) or None
      region: FR1/2/3/4 (str)
    Return:
      imgt numbered seq with gaps (str)
    '''
    # seq = seq.replace('-', '').replace('.', '')

    if region == "FR1":
        if len(seq) == 0: 
            return "." * 26
        if len(seq) == 25: 
            return seq[0:9] + "." + seq[9:25]
        if len(seq) == 26: 
            return seq
        return "X" * 26

    if region == "FR2":
        if len(seq) == 0:
            return "." * 17
        if len(seq) == 17:
            return seq
        return "X" * 17

    if region == "FR3":
        if len(seq) == 0:
            return "." * 39
        if len(seq) == 38:
            return seq[0:7] + "." + seq[7:38]
        if len(seq) == 39:
            return seq
        return "X" * 39

    if region == "FR4":
        if len(seq) == 0:
            return "." * 11
        if len(seq) == 11:
            return seq
        return "X" * 11   


def cdr_numbering(seq, region):
    '''
    IMGT unique numbering for CDRs
      CDR1          : 12  (5-12)
      CDR2          : 10  (0-10)
      CDR3 basic    : 13  (5-13)
      CDR3 extend 1 : 31 (14-31)
      CDR3 extend 2 : 51 (32-51)
      CDR3 extend 3 : 71 (52-71)
      CDR3 extend 4 : 91 (72-91)
    Arguments: 
      seq: aa sequence (str) or None
      region: CDR1/2 (str)
    Return:
      imgt numbered seq with gaps (str)
    '''
    # seq = seq.replace('-', '').replace('.', '')
    # seq_length = len(seq)

    # Gap orders in CDRs
    if region == 'CDR1':
        imgt_length = 12
        gaps = [7, 6, 8, 5, 9, 4, 10]

    if region == 'CDR2':
        imgt_length = 10
        gaps = [6, 5, 7, 4, 8, 3, 9, 2, 10, 1]

    if region == 'CDR3':
        imgt_length = 13
        gaps = [7, 8, 6, 9, 5, 10, 4, 11]
    if region == 'CDR3_extend_1':
        imgt_length = 31
        gaps = [16, 17, 15, 18, 14, 19, 13, 20, 12, 21, 11, 22, 10, 23, 9, 24, 8]
    if region == 'CDR3_extend_2':
        imgt_length = 51
        gaps = [26, 27, 25, 28, 24, 29, 23, 30, 22, 31, 21, 32, 20, 33, 19, 34, 18, 35, 17]
    if region == 'CDR3_extend_3':
        imgt_length = 71
        gaps = [36, 37, 35, 38, 34, 39, 33, 40, 32, 41, 31, 42, 30, 43, 29, 44, 28, 45, 27]
    if region == 'CDR3_extend_4':
        gaps = [46, 47, 45, 48, 44, 49, 43, 50, 42, 51, 41, 52, 40, 53, 39, 54, 38, 55, 37]
    
        
    # Add gaps to seq  
    if len(seq) == 0:
        imgt_seq = '.' * imgt_length
    elif len(seq) > imgt_length:
        # TODO:
        # May need to mark abnormal sequences that may arise from sequencing error.
        imgt_seq = seq
    else:
        # Find the positions to insert gaps
        n_gaps = imgt_length - int(len(seq))
        
        if n_gaps == 0: 
            imgt_seq =  seq
        else:
            gap_sites = gaps[:n_gaps]
            gaps_start = min(gap_sites) - 1
            # Insert gaps
            imgt_seq = seq[:gaps_start] + '.' * n_gaps + seq[gaps_start:]
    
    return imgt_seq

