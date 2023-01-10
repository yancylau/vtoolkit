# Determine V sequence type with v classifier


import pandas as pd
from dict import letter_dict
import pickle



# Load trained classifier
v_pkl = "models/v_classifier.pkl"
with open(v_pkl, 'rb') as file:
    v_classifier = pickle.load(file)


def determineType(imgt_fr2):
    '''
    Determine sequence type based on hallmarks on FR2
    Arguments: 
      imgt_fr2: imgt numbered FR2 (AA, str)
    Return:
      sequence type (str)
    '''
    if len(imgt_fr2) != 17: 
        return None
    
    hallmarks = pd.Series(imgt_fr2).str.split('', expand=True).iloc[:, [4,11,12,14]].apply(lambda x:x.map(letter_dict), axis=1)
    #hallmarks.columns = ["a42", "a49", "a50", "a52"] 
    predict = v_classifier.predict(hallmarks)[0]

    return "VHH" if predict == 1 else "VH"


# ## Test:
# imgt_fr2 = "MSWVRLSPGKGLEWVSA"
# imgt_fr2 = "MSWFRQAPGKQYEEVAR"
# imgt_fr2 = "MSWVRQAPGKEPEWIAW"
# imgt_fr2 = "MGWFRQAPGKEREAVAG"
# imgt_fr2 = "LAWFRQASGKEREGIAG"
# type = determineType(imgt_fr2)
# print(type)