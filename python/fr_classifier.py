# Fit Naive Bayes Clssifier to get FRs

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import CategoricalNB
from sklearn.metrics import accuracy_score, confusion_matrix, roc_curve, roc_auc_score, recall_score, precision_score
from sklearn.metrics import f1_score, classification_report
from sklearn import metrics
from dict import letter_dict
from dataset import fr1, fr2, fr3, fr4
import pickle



# Convert catergeorical variables to numbeical
fr1_encoded = fr1.apply(lambda x:x.map(letter_dict), axis=1) 
fr2_encoded = fr2.apply(lambda x:x.map(letter_dict), axis=1) 
fr3_encoded = fr3.apply(lambda x:x.map(letter_dict), axis=1)
fr4_encoded = fr4.apply(lambda x:x.map(letter_dict), axis=1)


# Select features and labels
fr1_x = fr1_encoded.loc[:,"a1":"a26"]
fr2_x = fr2_encoded.loc[:,"a39":"a55"]
fr3_x = fr3_encoded.loc[:,"a66":"a104"]
fr4_x = fr4_encoded.loc[:,"a118":"a128"]
fr1_y = fr1_encoded.loc[:,"type"]
fr2_y = fr2_encoded.loc[:,"type"]
fr3_y = fr3_encoded.loc[:,"type"]
fr4_y = fr4_encoded.loc[:,"type"]
# print(fr1_x)
# print(fr1_y)


# Train NB model for FR1, FR2 and FR3
def fit_nb_model(fr_x, fr_y):
    '''
    Train NB model for FR1, FR2, and FR3 seperately
    '''
    # Split dataset
    X = fr_x
    y = fr_y
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, random_state=123)
    
    # Fit model
    nb_classifier = CategoricalNB(alpha=1.0, fit_prior=True, class_prior=None)
    nb_classifier.fit(X_train, y_train)
    
    # Predict on test data
    predicts = nb_classifier.predict(X_test)
    
    accuracy = accuracy_score(y_test, predicts)
    print("accuracy:", accuracy)

    # ROC cruve
    roc_cruve = metrics.plot_roc_curve(nb_classifier, X_test, y_test)
    
    return nb_classifier


# Fit NB Clssifier models for FR1, FR2 and FR3
fr1_classifier = fit_nb_model(fr1_x, fr1_y)
fr2_classifier = fit_nb_model(fr2_x, fr2_y)
fr3_classifier = fit_nb_model(fr3_x, fr3_y)
fr4_classifier = fit_nb_model(fr4_x, fr4_y)


# Save trained models to file
fr1_pkl = "../models/fr1.pkl"
fr2_pkl = "../models/fr2.pkl"
fr3_pkl = "../models/fr3.pkl"
fr4_pkl = "../models/fr4.pkl"

with open(fr1_pkl, 'wb') as file:
    pickle.dump(fr1_classifier, file)
with open(fr2_pkl, 'wb') as file:
    pickle.dump(fr2_classifier, file)
with open(fr3_pkl, 'wb') as file:
    pickle.dump(fr3_classifier, file)
with open(fr4_pkl, 'wb') as file:
    pickle.dump(fr4_classifier, file)
