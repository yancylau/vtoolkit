# Extract variable (V) region from query sequence


from numbering import fr_numbering, cdr_numbering
from get_fr import getFR, getTruncatedFR



def getV(seq):
    '''
    Get FRs and CDRs from query sequence
    Arguments:
      seq: aa sequence (str)
    Returns:
      vregion: seqs and borders of FRs and CDRs (dict)
    '''
    vregion = {"fr1_start": None, 
               "fr1_end": None, 
               "fr1": "",
               "fr1_imgt": "",
               "fr2_start": None, 
               "fr2_end": None, 
               "fr2": "",
               "fr2_imgt": "",
               "fr3_start": None, 
               "fr3_end": None, 
               "fr3": "",
               "fr3_imgt": "",
               "fr4_start": None,
               "fr4_end": None,
               "fr4": "",
               "fr4_imgt": "",
               "cdr1_start": None, 
               "cdr1_end": None, 
               "cdr1": "",
               "cdr1_imgt": "",
               "cdr2_start": None, 
               "cdr2_end": None, 
               "cdr2": "",
               "cdr2_imgt": "",
               "cdr3_start": None, 
               "cdr3_end": None, 
               "cdr3": "",
               "cdr3_imgt": "",
               "v_aa": "", 
               "v_imgt_aa": ""}

    # TODO: This criterion is arbitrary. FWD-aa REV-aa need to be optimized.
    if len(seq) < 40: 
        return vregion

    # Detect FRs in the query sequence
    fr1 = getFR(seq, "FR1")
    fr2 = getFR(seq, "FR2")
    fr3 = getFR(seq, "FR3")
    fr4 = getFR(seq, "FR4")

    # If FR not found, further search subseqs at heads and tails
    if len(fr1["seq"]) == 0:
        fr1 = getTruncatedFR(seq, "FR1", "tail")
    if len(fr2["seq"]) == 0:
        fr2 = getTruncatedFR(seq, "FR2", "tail")
    if len(fr3["seq"]) == 0:
        fr3 = getTruncatedFR(seq, "FR3", "tail")
    if len(fr4["seq"]) == 0:
        fr4 = getTruncatedFR(seq, "FR4", "tail")

    if (fr1["seq"] == fr2["seq"] == fr3["seq"] == fr4["seq"] == ""): 
        return vregion

    vregion['fr1_start'] = fr1['start'] 
    vregion['fr1_end'] = fr1['end']
    vregion['fr1'] = fr1['seq']
    vregion['fr2_start'] = fr2['start']
    vregion['fr2_end'] = fr2['end']
    vregion['fr2'] = fr2['seq']
    vregion['fr3_start'] = fr3['start']
    vregion['fr3_end'] = fr3['end']
    vregion['fr3'] = fr3['seq']
    vregion['fr4_start'] = fr4['start']
    vregion['fr4_end'] = fr4['end']
    vregion['fr4'] = fr4['seq']
    
    # Get CDRs in the query sequence
    if fr1['end'] and fr2['start']:
        vregion['cdr1_start'] = fr1['end'] + 1
        vregion['cdr1_end'] = fr2['start'] - 1
        vregion['cdr1'] = seq[vregion['cdr1_start']-1:vregion['cdr1_end']]
    if fr2['end'] and fr3['start']:
        vregion['cdr2_start'] = fr2['end'] + 1
        vregion['cdr2_end']= fr3['start'] - 1
        vregion['cdr2'] = seq[vregion['cdr2_start']-1:vregion['cdr2_end']]
    if fr3['end'] and fr4['start']:
        vregion['cdr3_start'] = fr3['end'] + 1
        vregion['cdr3_end']= fr4['start'] - 1
        vregion['cdr3'] = seq[vregion['cdr3_start']-1:vregion['cdr3_end']]
    # AA sequence of V region
    vregion['v_aa'] =  vregion['fr1'] + vregion['cdr1'] + vregion['fr2'] + vregion['cdr2'] + vregion['fr3'] + vregion['cdr3'] + vregion['fr4'] 

    # Add IMGT numbered FR1/FR2/FR3 and CDR1/2
    vregion['fr1_imgt'] = fr_numbering(vregion['fr1'], "FR1")
    vregion['fr2_imgt'] = fr_numbering(vregion['fr2'], "FR2")
    vregion['fr3_imgt'] = fr_numbering(vregion['fr3'], "FR3")
    vregion['fr4_imgt'] = fr_numbering(vregion['fr4'], "FR4")
    vregion['cdr1_imgt'] = cdr_numbering(vregion['cdr1'], "CDR1")
    vregion['cdr2_imgt'] = cdr_numbering(vregion['cdr2'], "CDR2")

    # Add IMGT numbered CDR3
    if len(vregion['cdr3']) <= 13:
        vregion["cdr3_length"] = "basic"
        vregion['cdr3_imgt'] = cdr_numbering(vregion['cdr3'], "CDR3")
    elif len(vregion['cdr3']) <= 31:
        vregion["cdr3_length"] = "extend_1"
        vregion['cdr3_imgt'] = cdr_numbering(vregion['cdr3'], "CDR3_extend_1")
    elif len(vregion['cdr3']) <= 51:
        vregion["cdr3_length"] = "extend_2"
        vregion['cdr3_imgt'] = cdr_numbering(vregion['cdr3'], "CDR3_extend_2")
    elif len(vregion['cdr3']) <= 71:
        vregion["cdr3_length"] = "extend_3"
        vregion['cdr3_imgt'] = cdr_numbering(vregion['cdr3'], "CDR3_extend_3")
    elif len(vregion['cdr3']) <= 91:
        vregion["cdr3_length"] = "extend_4"
        vregion['cdr3_imgt'] = cdr_numbering(vregion['cdr3'], "CDR3_extend_4")
    else:
        vregion["cdr3_length"] = "other"
        vregion['cdr3_imgt'] = "." * 13

    # Add IMGT numbered V region
    vregion['v_imgt_aa'] =  vregion['fr1_imgt'] + vregion['cdr1_imgt'] + vregion['fr2_imgt'] + vregion['cdr2_imgt'] + vregion['fr3_imgt'] + vregion['cdr3_imgt'] + vregion['fr4_imgt'] 

    
    return vregion



# # Test:
# seq = "QVQLVESGGGLVQPGGSLRLSCAASGFTLDYYAIGWFRQAPGKEREGVSCISSSDGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA"
# vregion = getV(seq)
# print(vregion)