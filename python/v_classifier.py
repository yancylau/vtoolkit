# Fit V classifier with Naive Bayes Clssifier
# Let 4 classical hallmarks on FR2 as features



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


# Select features and labels
frs_x = frs_encoded.loc[:, ("a42", "a49", "a50", "a52")]
frs_y = frs_encoded.loc[:, "type"]


# Function: Fit model with NB classifer
def fit_nb_model(fr_x, fr_y):
    '''
    Fit NB model for V classifier
    '''
    # Split dataset into training and test subsets
    # To avoid waning, use "fr_x.values" instead of "fr_x"
    # Ref: https://stackoverflow.com/questions/69326639/sklearn-warning-valid-feature-names-in-version-1-0
    X = fr_x
    y = fr_y
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=123)
    
    # Fit model
    nb_classifier = CategoricalNB(alpha=1.0, fit_prior=True, class_prior=None)
    nb_classifier.fit(X_train, y_train)
    
    # Predict on test data
    predicts = nb_classifier.predict(X_test)

    # Performance
    accuracy = accuracy_score(y_test, predicts)
    
    return nb_classifier


# Fit NB Clssifier model
v_classifier = fit_nb_model(frs_x.values, frs_y)


# Save fitted models to file
v_pkl = "models/v_classifier.pkl"
with open(v_pkl, 'wb') as file:
    pickle.dump(v_classifier, file)
