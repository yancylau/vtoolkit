# Gererate training and test dataset for FR1/2/3/4 detectors and V classifier


import pandas as pd
import numpy as np
import string
from dict import aas



## 0. Read into IGHV sequence dataset
v_seqs = pd.read_csv('data/v_reference_extend.csv')


## 1. Gererate dataset for FR1/2/3 detector
# Positive dataset
fr1_positive = v_seqs["FR1"].str.split('', expand=True).iloc[:, 1:27].replace({"X": None, ".": None, "-": None })
fr2_positive = v_seqs["FR2"].str.split('', expand=True).iloc[:, 1:18].replace({"X": None, ".": None, "-": None })
fr3_positive = v_seqs["FR3"].str.split('', expand=True).iloc[:, 1:40].replace({"X": None, ".": None, "-": None })
fr1_positive[27] = "positive"
fr2_positive[18] = "positive"
fr3_positive[40] = "positive"

# Negative dataset
# TODO: non-IGHV sequences from refseq
# letters = list(string.ascii_uppercase)
# aas = letters
fr1_negative = pd.DataFrame(np.random.choice(aas, size=(10000, 26)), columns=list(range(1,27)))
fr2_negative = pd.DataFrame(np.random.choice(aas, size=(10000, 17)), columns=list(range(1,18)))
fr3_negative = pd.DataFrame(np.random.choice(aas, size=(10000, 39)), columns=list(range(1,40)))
fr1_negative[27] = "negative"
fr2_negative[18] = "negative"
fr3_negative[40] = "negative"

# Concatenate positive dataset with negative dataset
# NOTE: Remove the null site in FR1 and FR3
fr1 = pd.concat([fr1_positive, fr1_negative], axis=0)
fr2 = pd.concat([fr2_positive, fr2_negative], axis=0)
fr3 = pd.concat([fr3_positive, fr3_negative], axis=0)

# Set columan names: features and label
fr1.columns =['a{}'.format(i) for i in range(1,27)] + ["type"]
fr2.columns =['a{}'.format(i) for i in range(39,56)] + ["type"]
fr3.columns = ['a{}'.format(i) for i in range(66,105)] + ["type"]

# Remove columns with gaps
fr1 = fr1.drop(["a10"], axis=1)
fr3 = fr3.drop(["a73"], axis=1)


## 2. Gererate dataset for FR4 classifier
ighj = pd.read_table("data/fr4.tsv", header = None)
ighj.columns =  ['id', 'ighj', 'fr4']
ighj = ighj[ighj['fr4'].str.len() == 11]

# Positive
fr4_positive =  ighj['fr4'].str.split('', expand=True).iloc[:, 1:12].replace({"X": None, ".": None, "-": None, np.nan: None})
fr4_positive[13] = "positive"

# Negative dataset
letters = list(string.ascii_uppercase)
fr4_negative = pd.DataFrame(np.random.choice(list(string.ascii_uppercase), size=(10000, 11)), columns=list(range(1,12)))
fr4_negative[13] = "negative"

# Merge positive with negative and set columan names
fr4 = pd.concat([fr4_positive, fr4_negative], axis=0)
fr4.columns =['a{}'.format(i) for i in range(118,129)] + ["type"]


## 3. Gererate trainng and test dataset for VH/VHH classifier
frs = v_seqs["Full"].str.split('', expand=True).iloc[:, np.r_[1:27, 39:56, 66:105]].replace({"X": None, ".": None, "-": None })
frs.columns =  ['a{}'.format(i) for i in list(range(1,27)) + list(range(39,56)) + list(range(66,105))]
frs['type'] = v_seqs["Type"]




# # Test
# print("FR1", fr1)
# print("FR2", fr2.columns)
# print("FR3", fr3.columns)
# print("FR4", fr4.columns)