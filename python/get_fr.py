# Use FR1/2/3/4 detectors to extract target sequence from inquery


import pandas as pd
import numpy as np
from functions import get_heads, get_tails
from dict import letter_dict
import pickle



# Load trained models
fr1_pkl = "models/fr1.pkl"
fr2_pkl = "models/fr2.pkl"
fr3_pkl = "models/fr3.pkl"
fr4_pkl = "models/fr4.pkl"
with open(fr1_pkl, 'rb') as file:
    fr1_detector = pickle.load(file)
with open(fr2_pkl, 'rb') as file:
    fr2_detector = pickle.load(file)
with open(fr3_pkl, 'rb') as file:
    fr3_detector = pickle.load(file)
with open(fr4_pkl, 'rb') as file:
    fr4_detector = pickle.load(file)


def getFR(seq, region):
    '''
    To get FR seq from query seq
    Arguments: 
      seq: aa sequence (str)
    Return:
      result: FR region (dict)
    '''
    fr = {"start": None, 
          "end": None, 
          "seq": "",
          "seq_extend": ""}

    # Assign FR arguments
    if region == "FR1":
        region_length = 25
        feature_names = pd.Series(['a{}'.format(i) for i in range(1,27)]).drop(9)
        fr_detector = fr1_detector
    if region == "FR2":
        region_length = 17
        feature_names = pd.Series(['a{}'.format(i) for i in range(39,56)])
        fr_detector = fr2_detector
    if region == "FR3":
        region_length = 38
        feature_names = pd.Series(['a{}'.format(i) for i in range(66,105)]).drop(7)
        fr_detector = fr3_detector
    if region == "FR4":
        region_length = 11
        feature_names = pd.Series(['a{}'.format(i) for i in range(118,129)])
        fr_detector = fr4_detector

    # If query sequence length is shorter than FR length, return None
    if len(seq) < region_length: 
        return fr

    # Determine FR and its border from subseqs
    positives = []
    probs = []
    predict = 0
    n_subseqs = len(seq) - region_length + 1
    for i in range(0, n_subseqs):
        subseq = seq[i:i+region_length]
        subseq = [letter_dict.get(aa) for aa in subseq] 
        subseq = np.array(subseq).reshape(1, -1)
        predict = fr_detector.predict(subseq)[0]
        if predict == 1:
            positives.append(i)
            probs.append(fr_detector.predict_proba(subseq)[0,1])
    if len(positives) > 0:
        start = positives[probs.index(max(probs))] + 1
        end = positives[probs.index(max(probs))] + region_length
        fr['start'] = start
        fr['end'] = end
        fr['seq'] = seq[start-1:end]

    return fr


def getTruncatedFR(seq, region, end = "head"):
    '''
    To get truncated FR at the end (head/tail) of the seq
    Arguments: 
      seq: aa sequence (str)
      region: region name (str)
      end: head/end (str)
    Return:
      fr: FR subregion (dict)
    '''
    truncated_fr = {'start': None, 
                    'end': None, 
                    'seq': "", 
                    'seq_extend': ""}

    # Assign arguments for each FR
    if region == "FR1":
        region_length = 25
        feature_names = pd.Series(['a{}'.format(i) for i in range(1,27)]).drop(9)
        fr_detector = fr1_detector
    if region == "FR2":
        region_length = 17
        feature_names = pd.Series(['a{}'.format(i) for i in range(39,56)])
        fr_detector = fr2_detector
    if region == "FR3":
        region_length = 38
        feature_names = pd.Series(['a{}'.format(i) for i in range(66,105)]).drop(7)
        fr_detector = fr3_detector
    if region == "FR4":
        region_length = 11
        feature_names = pd.Series(['a{}'.format(i) for i in range(118,129)])
        fr_detector = fr4_detector
        
    if end == "head": 
        subseqs = get_heads(seq, region_length)
        #subseqs.columns = feature_names
        encoded_subseqs = subseqs.replace({"X": None, ".": None, "-": None }).apply(lambda x:x.map(letter_dict), axis=1) 
        predicts = fr_detector.predict(encoded_subseqs)
        positives = np.where(predicts==1)[0]

        if positives.size == 0: 
            # If FR is NOT found, return None
            return truncated_fr
        else:
            proba = fr_detector.predict_proba(encoded_subseqs)
            max_proba_index = np.argmax(proba[:,1])
            truncated_fr['start'] = 1
            truncated_fr['end'] = region_length - max_proba_index - 1
            truncated_fr['seq'] = seq[:region_length - max_proba_index - 1]
            truncated_fr['seq_extend'] = subseqs.iloc[max_proba_index].str.cat(sep="")


    if end == "tail": 
        subseqs = get_tails(seq, region_length)
        #subseqs.columns = feature_names
        encoded_subseqs = subseqs.replace({"X": None, ".": None, "-": None }).apply(lambda x:x.map(letter_dict), axis=1) 
        predicts = fr_detector.predict(encoded_subseqs)
        positives = np.where(predicts==1)[0]

        if positives.size == 0: 
            # If FR is NOT found, return None
            return truncated_fr
        else:
            proba = fr_detector.predict_proba(encoded_subseqs)
            max_proba_index = np.argmax(proba[:,1])
            truncated_fr['start'] = len(seq) - region_length +  max_proba_index + 2
            truncated_fr['end'] = len(seq)
            truncated_fr['seq'] = subseqs.iloc[max_proba_index].str.cat(sep="")
            # truncated_fr['seq'] = seq[len(seq) - region_length +  max_proba_index + 2:len(seq)]
            # truncated_fr['seq_truncated'] = subseqs.iloc[max_proba_index].str.cat(sep="")
            
    return(truncated_fr)




# ### Test:
# seq = "AAEVQLVESGGGSVQVGGSLNLSCPVSGYTRGHYLMAWFRQVAGKEREGVAGIDSENSAVYAGSVKGRFTISLDNAKNTLYLRMTSLKPEDTGMYYCAADRGLLFRPLYSDTYKNWGQGTQVTVSSKKKK"
# fr1 =  getFR(seq, "FR1")
# fr2 =  getFR(seq, "FR2")
# fr3 =  getFR(seq, "FR3")
# fr4 =  getFR(seq, "FR4")
# # truncated_fr4 =  getTruncatedFR(seq, "FR4", end = "tail")
# print("FR1", fr1, len(fr1["seq"]))
# print("FR2", fr2)
# print("FR3", fr3)
# print("FR4", fr4)
# # # print(truncated_fr4)



