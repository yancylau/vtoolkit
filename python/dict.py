# Encode categorical features. 
# Dictionary to convert alphabet letters to numbers


import string



# Alphabet letters and labels
letters = list(string.ascii_uppercase)

keys = [None] + letters + ["negative", "positive", "VH", "VHH"]
values = [0] + list(range(1,27)) + [0,1,0,1]

letter_dict = dict(zip(keys, values))


# Amino acid list, without ambigous amino acids (BJOUXZ)
aas = ['A', 
       'C', 
       'D', 
       'E', 
       'F', 
       'G', 
       'H', 
       'I', 
       'K', 
       'L', 
       'M', 
       'N', 
       'P', 
       'Q', 
       'R', 
       'S', 
       'T', 
       'V', 
       'W', 
       'Y']




# # Exclude ambious amino acids: BJOUXZ
# letter_dict = {None: 0, 
#                'A': 1, 
#                'C': 2, 
#                'D': 3, 
#                'E': 4, 
#                'F': 5, 
#                'G': 6, 
#                'H': 7, 
#                'I': 8, 
#                'K': 9, 
#                'L': 10, 
#                'M': 11, 
#                'N': 12, 
#                'P': 13, 
#                'Q': 14, 
#                'R': 15, 
#                'S': 16, 
#                'T': 17, 
#                'V': 18, 
#                'W': 19, 
#                'Y': 20, 
#                'negative': 0, 
#                'positive': 1,
#                'VH': 0,
#                'VHH': 1}