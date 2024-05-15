from MVPA_func import SpherePoints, SphereValue, LoadAllNii, MVPA_FeatureExtraction, MVPA_Classification, LoopiFeature, multi_run_wrapper

import numpy as np
import pandas as pd
import os.path as op
import matplotlib.pyplot as plt
import nibabel as nib
import itertools

from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
#from sklearn.svm import LinearSVC

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.pipeline import Pipeline
from sklearn.multioutput import MultiOutputClassifier

def multi_run_wrapper(args):
    return LoopiFeature(*args)