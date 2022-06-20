# Use FR detectors to extract subregions from aa sequence

from numbering import fr_numbering, cdr_numbering
from get_fr import getFR, getExtendFR




def getSubregions(seq):
    '''
    Get FRs and CDRs from query sequence
    Arguments:
      seq: aa sequence (str)
    Returns:
      subregion: seqs and borders of FRs and CDRs (dict)
    '''
    subregion = {"fr1_start": None, 
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
        return subregion

    # Detect FRs in the query sequence
    fr1 = getFR(seq, "FR1")
    fr2 = getFR(seq, "FR2")
    fr3 = getFR(seq, "FR3")
    fr4 = getFR(seq, "FR4")

    # If FR not found, further search subseqs at heads and tails
    if len(fr1["seq"]) == 0:
        fr1 = getExtendFR(seq, "FR1", "tail")
    if len(fr2["seq"]) == 0:
        fr2 = getExtendFR(seq, "FR2", "tail")
    if len(fr3["seq"]) == 0:
        fr3 = getExtendFR(seq, "FR3", "tail")
    if len(fr4["seq"]) == 0:
        fr4 = getExtendFR(seq, "FR4", "tail")

    if (fr1["seq"] == fr2["seq"] == fr3["seq"] == fr4["seq"] == ""): 
        return subregion

    subregion['fr1_start'] = fr1['start'] 
    subregion['fr1_end'] = fr1['end']
    subregion['fr1'] = fr1['seq']
    subregion['fr2_start'] = fr2['start']
    subregion['fr2_end'] = fr2['end']
    subregion['fr2'] = fr2['seq']
    subregion['fr3_start'] = fr3['start']
    subregion['fr3_end'] = fr3['end']
    subregion['fr3'] = fr3['seq']
    subregion['fr4_start'] = fr4['start']
    subregion['fr4_end'] = fr4['end']
    subregion['fr4'] = fr4['seq']
    
    # Determine CDRs in the query sequence
    if fr1['end'] and fr2['start']:
        subregion['cdr1_start'] = fr1['end'] + 1
        subregion['cdr1_end'] = fr2['start'] - 1
        subregion['cdr1'] = seq[subregion['cdr1_start']-1:subregion['cdr1_end']]
    if fr2['end'] and fr3['start']:
        subregion['cdr2_start'] = fr2['end'] + 1
        subregion['cdr2_end']= fr3['start'] - 1
        subregion['cdr2'] = seq[subregion['cdr2_start']-1:subregion['cdr2_end']]
    if fr3['end'] and fr4['start']:
        subregion['cdr3_start'] = fr3['end'] + 1
        subregion['cdr3_end']= fr4['start'] - 1
        subregion['cdr3'] = seq[subregion['cdr3_start']-1:subregion['cdr3_end']]
    # AA sequence of V region
    subregion['v_aa'] =  subregion['fr1'] + subregion['cdr1'] + subregion['fr2'] + subregion['cdr2'] + subregion['fr3'] + subregion['cdr3'] + subregion['fr4'] 

    # Add IMGT numbered sequences
    subregion['fr1_imgt'] = fr_numbering(subregion['fr1'], "FR1")
    subregion['fr2_imgt'] = fr_numbering(subregion['fr2'], "FR2")
    subregion['fr3_imgt'] = fr_numbering(subregion['fr3'], "FR3")
    subregion['fr4_imgt'] = fr_numbering(subregion['fr4'], "FR4")
    subregion['cdr1_imgt'] = cdr_numbering(subregion['cdr1'], "CDR1")
    subregion['cdr2_imgt'] = cdr_numbering(subregion['cdr2'], "CDR2")
    subregion['cdr3_imgt'] = cdr_numbering(subregion['cdr3'], "CDR3")
    subregion['v_imgt_aa'] =  subregion['fr1_imgt'] + subregion['cdr1_imgt'] + subregion['fr2_imgt'] + subregion['cdr2_imgt'] + subregion['fr3_imgt'] + subregion['cdr3_imgt'] + subregion['fr4_imgt'] 

    # Add sequence type (VH/VHH)
    # if subregion["imgt_fr2"]:
    #     subregion["type"] = determine_type(subregion["imgt_fr2"])
    # else:
    #     subregion['type'] = None
    
    return subregion




### Test:
# seq = "MAHVQLVESGGGSVQAGGSLRLSCVASGNTKCMAWFRQAPGKEREGVATIHHRSLGALTYYADSVKGQFTISRDYAKNTLYLQLNSLKTEDTAMYYCAKDSPRGQWTLIATILDEYNYWGQGTQVTVS"
# seq = "MADVQLVESGGGLVQPGGSLRLSCAAFGFTFSSYDMSWVRRAPGKGLEWVSAINSVGSSTYYTDSVNGRFTISRDNAKNTLYLQMNSLKTEDTAVYYCATGTSGRWYDRDSGYWGQGTQVTVS"
# subregion = getSubregions(seq)
# print(subregion)