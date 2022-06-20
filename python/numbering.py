# IMGT unique numbering for FRs and CDRs


def fr_numbering(seq, region):
    '''
    IMGT unique numbering for FRs
    FR1:   1 -  26
    FR2:  39 -  55
    FR3:  66 - 104
    FR4: 118 - 128
    Usually, 
      FR1 is 25 aa with one gap in position 10;
      FR2 is 12 aa withnot gaps;
      FR3 is 38 aa with one gap in position 8; and
      FR4 is 11 aa without gaps

    Arguments: 
      seq: aa sequence (str) or None
      region: FR1/2/3/4 (str)
    Return:
      imgt numbered seq with gaps (str)
    '''
    # seq = seq.replace('-', '').replace('.', '')
    # seq_length = len(seq)

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
    CDR1_length: 5-12
    CDR2_length: 0-10
    CDR3_length: 5-13
    CDR3_length (extend 1): 14-31
    CDR3_length (extend 2): 32-51

    Ref: http://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html#table1

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
        
    # Add gaps to seq  
    if len(seq) == 0:
        imgt_seq = '.' * imgt_length
    elif len(seq) > imgt_length:
        # TODO: May need to mark abnormal sequences 
        # In fact, this sequence may arised by sequencing error
        imgt_seq = seq
    else:
        # Find the positions to insert gaps
        n_gaps = imgt_length - int(len(seq))

        if n_gaps == 0: 
            imgt_seq =  seq
        else:
            gap_sites = gaps[:n_gaps]
            gaps_start = min(gap_sites) - 1
            # Insert gaps to seq
            imgt_seq = seq[:gaps_start] + '.' * n_gaps + seq[gaps_start:]
    
    return imgt_seq

