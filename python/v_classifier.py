# Fit Naive Bayes Clssifier to distinguish/seperate VHHs from VHs

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import CategoricalNB
from sklearn.metrics import accuracy_score, confusion_matrix, roc_curve, roc_auc_score, recall_score, precision_score
from sklearn.metrics import f1_score, classification_report
from dict import letter_dict
from dataset import frs
import pickle

# Convert catergeorical variables to numbeical
frs_encoded = frs.apply(lambda x:x.map(letter_dict), axis=1) 



# Fit model with hallmarks only


# Select features and labels
# frs_x = frs_encoded.loc[:,("a42", "a49","a50", "a52")]
frs_x = frs_encoded.loc[:,("a42", "a49","a50", "a52")]
frs_y1 = frs_encoded.loc[:,"type"]
frs_y2 = frs_encoded.loc[:,"type_adjusted"]




# print(fr1_x)
# print(fr1_y)


# Function: Train NB model
def fit_nb_model(fr_x, fr_y):
    '''
    Train NB model for FR1, FR2, and FR3 seperately
    '''

    # Split dataset into train and test subsets
    X = fr_x
    y = fr_y
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)
    
    # Fit model
    nb_classifier = CategoricalNB(alpha=1.0, fit_prior=True, class_prior=None)
    nb_classifier.fit(X_train, y_train)
    
    # Predict on test data
    predicts = nb_classifier.predict(X_test)
    accuracy = accuracy_score(y_test, predicts)
    #print("accuracy:", accuracy)
    
    return nb_classifier


# Fit NB Clssifier model
# v_classifier_hallmarks = fit_nb_model(frs_x, frs_y1)
v_classifier_hallmarks = fit_nb_model(frs_x, frs_y2)

# Return the final classifier
v_classifier = v_classifier_hallmarks


# Save trained models to file
v_pkl = "../models/v_classifier.pkl"
with open(v_pkl, 'wb') as file:
    pickle.dump(v_classifier, file)
