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

def SpherePoints(center:  np.array, radius: int, voxel_dims: tuple) -> np.array :
    '''This function generates sphere around around a center and predefined sphere'''
    '''It avoids incomplete spheres'''
    '''It returns the coordinates of the points aroud r distance from center '''

    
    #outer bounding box
    x0 = center[0]-radius
    x1 = center[0]+radius
    y0 = center[1]-radius
    y1 = center[1]+radius
    z0 = center[2]-radius
    z1 = center[2]+radius
    cnt=1
    sphere=np.array([])
    spherex=[]
    spherey=[]
    spherez=[]
    if ((x0>=0) & (y0>=0) & (z0>=0) &(x1<voxel_dims[0]) & (y1<voxel_dims[1]) & (z1<voxel_dims[2])):
    #     print('X:',x0,x1)
    #     print('Y:',y0,y1)
    #     print('Z:',z0,z1)
        #add 1 to upper bounds to range is inclusive

        for coords in itertools.product(range(x0,x1+1), range(y0,y1+1), range(z0,z1+1)):
            if radius**2 >= sum((center - coords)**2): #euclidean distance is smaller than radius
                spherex.append(coords[0])
                spherey.append(coords[1])
                spherez.append(coords[2])
#                 if cnt==1: 
#                     sphere=np.array([coords])
#                     cnt=0
#                 else:
#                     sphere=np.append(sphere,[coords],axis=0)
  #  print(sphere)
#    print(spherex,spherey,spherez)
    sphere=np.array([[spherex],[spherey],[spherez]])
    return sphere
    
def IncompleteSpherePoints(center:  np.array, radius: int, voxel_dims: tuple) -> np.array :
    '''This function generates sphere around around a center and predefined sphere'''
    '''It DOES NOT avoid incomplete spheres'''
    '''It returns the coordinates of the points aroud r distance from center '''

    
    #outer bounding box
    x0 = center[0]-radius
    x1 = center[0]+radius
    y0 = center[1]-radius
    y1 = center[1]+radius
    z0 = center[2]-radius
    z1 = center[2]+radius
    cnt=1
    sphere=np.array([])
    spherex=[]
    spherey=[]
    spherez=[]
#    if ((x0>=0) & (y0>=0) & (z0>=0) &(x1<voxel_dims[0]) & (y1<voxel_dims[1]) & (z1<voxel_dims[2])):
    #     print('X:',x0,x1)
    #     print('Y:',y0,y1)
    #     print('Z:',z0,z1)
        #add 1 to upper bounds to range is inclusive

    for coords in itertools.product(range(x0,x1+1), range(y0,y1+1), range(z0,z1+1)):
        if radius**2 >= sum((center - coords)**2): #euclidean distance is smaller than radius
            spherex.append(coords[0])
            spherey.append(coords[1])
            spherez.append(coords[2])
#                 if cnt==1: 
#                     sphere=np.array([coords])
#                     cnt=0
#                 else:
#                     sphere=np.append(sphere,[coords],axis=0)
  #  print(sphere)
#    print(spherex,spherey,spherez)
    sphere=np.array([[spherex],[spherey],[spherez]])
    return sphere
    
def SphereValue(myPoints:  np.array, data: np.array) -> np.array :
    '''Extract the value of the sphere points from data '''
    '''direct indexing is faster, so I'm not using this code'''
    sphereVal=np.array([])
    if myPoints.shape!=(0,):
        cnt=1
        for coord0,coord1,coord2 in myPoints:
            if cnt==1: 
                sphereVal=np.array(data[coord0,coord1,coord2])
                cnt=0
            else:
                sphereVal=np.append(sphereVal,data[coord0,coord1,coord2])   

    return sphereVal


def LoadAllNii(ZstatInfoSubj,CLconditions):
    '''Load Nii data and keep the value in Alldata'''
    cntAllDataFirst=1
    for iRun in ZstatInfoSubj['Run']:
        MyZstat=ZstatInfoSubj.loc[(ZstatInfoSubj['Run']==iRun)]
        MyZstat.reset_index(inplace = True, drop = True)
        for iCond,CondName in enumerate(CLconditions):
             for MyInd,MyNii in enumerate(MyZstat[CondName][0]):
                MyFile=MyZstat['Folder'][0]+MyNii
                data = nib.load(MyFile).get_data()
                if cntAllDataFirst==1:
                    Alldata=np.array([data])
                    cntAllDataFirst=0
                else:
                    Alldata=np.append(Alldata,np.array([data]), axis=0)
    return Alldata
    
def MVPA_FeatureExtraction(ZstatInfoSubj,CLconditions,Alldata, myPoints):
    '''In this function features for a classifer is extracted based ont he myPoints provided'''
    
    cntAllData=0
    X=np.empty((0, myPoints.shape[-1]), int) # the number of voxels in the spot light search
    Y=np.array([])    
    Gr=np.array([]) 
    for iRun in ZstatInfoSubj['Run']:
        #print(iRun)
        MyZstat=ZstatInfoSubj.loc[(ZstatInfoSubj['Run']==iRun)]
        MyZstat.reset_index(inplace = True, drop = True)
        #MyZstat=MyZstat.loc[0]
        X_run=np.empty((0, myPoints.shape[-1]), int)
        Y_run=np.array([])
        Gr_run=np.array([])
        for iCond,CondName in enumerate(CLconditions):
            #print(MyZstat[CondName][0])
            X_tmp = np.zeros((len(MyZstat[CondName][0]), myPoints.shape[-1]))
            Y_tmp = np.zeros(len(MyZstat[CondName][0]))
            for MyInd,MyNii in enumerate(MyZstat[CondName][0]):
                MyFile=MyZstat['Folder'][0]+MyNii
#                     data = nib.load(MyFile).get_data()
#                     Alldata=np.array([data])
#                     Alldata=np.append(Alldata,np.array([data]), axis=0)
                data=Alldata[cntAllData,:,:,:]
                cntAllData+=1

                #sphereVal = data[myPoints[0],myPoints[1],myPoints[2]] 
                sphereVal = data[np.int_(myPoints[0]),np.int_(myPoints[1]),np.int_(myPoints[2])] 

                #data = data.ravel()
                X_tmp[MyInd,:] = sphereVal
                Y_tmp[MyInd]=int(iCond)+1

            X_run=np.concatenate((X_run,X_tmp),axis = 0)
            Y_run=np.concatenate((Y_run,Y_tmp),axis = 0)
            Gr_run=np.concatenate((Gr_run,np.repeat(iRun,Y_tmp.shape)),axis = 0)
            #print(Y.shape)
        X=np.concatenate((X,X_run),axis = 0)
        Y=np.concatenate((Y,Y_run),axis = 0)
        Gr=np.concatenate((Gr,Gr_run),axis = 0)
#     print(cntAllData)
    return X, Y, Gr
    
    
    
def MVPA_Classification(X,Y,Gr, pipeSpotLight):
    '''Performe classification on voxed data'''
    SpotPerformance = 0;
    NRun=int(max(Gr))
    if (np.any(X)): # perform classifier if values of X are not all zero
        # print(Gr)
        gkf = GroupKFold(n_splits=NRun) # leave one out, initialize GroupKFold with Nrun splits
        performance_this_fold = np.zeros(gkf.n_splits)
        for i_fold, (train_idx, test_idx) in enumerate(gkf.split(X=X, y=Y, groups=Gr)):
            # print("Indices of our test-samples: %r" % test_idx.tolist())
            # print("... which correspond to following runs: %r" % Gr[test_idx].tolist(), '\n')
            # Implement your ToDo here!
            X_test, X_train=X[test_idx], X[train_idx]
            Y_test, Y_train=Y[test_idx], Y[train_idx]

            pipeSpotLight.fit(X_train, Y_train)
            preds = pipeSpotLight.predict(X_test)
    #                 print('TRUE:',Y_test)
    #                 print('Pred:',preds)
    #                 print("Accuracy test: %.2f" % (preds == Y_test).mean())

    #        performance = roc_auc_score(Y_test,  pipe.predict_proba(X_test), multi_class='ovr') # 'ovr'(One-vs-rest): Computes the AUC of each class against the rest, This treats the multiclass case in the same way as the multilabel case. Sensitive to class imbalance even when average == 'macro', because class imbalance affects the composition of each of the ‘rest’ groupings.
            # in the pairwise AUC calculation, if a sample is labeled wrongly as a another class (not in paired), it will be dropped from the calculation, resulting in a higher rate
            performance = roc_auc_score(Y_test,  pipeSpotLight.predict_proba(X_test),  multi_class='ovo', average = 'macro') # 'ovo' (One-vs-one) Computes the average AUC of all possible pairwise combinations of classes. Insensitive to class imbalance when average == 'macro'.

            performance_this_fold[i_fold] = performance

        SpotPerformance=np.mean(performance_this_fold)
        
    return SpotPerformance
    
    
def LoopiFeature(iFeature: int, AllPoints: np.ndarray,ZstatInfoSubj: pd.core.frame.DataFrame ,CLconditions: list, Alldata: np.ndarray, pipeSpotLight: Pipeline):
    myPoints = AllPoints[iFeature] 
    # myPoints = SpherePoints(np.array([30,30,15]) , radius) # sample good point

    [X, Y, Gr]=MVPA_FeatureExtraction(ZstatInfoSubj,CLconditions,Alldata, myPoints) # extract data for classifer
    SpotPerformance=MVPA_Classification(X,Y,Gr, pipeSpotLight)
    return SpotPerformance
    
def multi_run_wrapper(args):
    return LoopiFeature(*args)
    
def LoopiFeatureDict(iFeature: int, DataArgDict: dict):
    '''This function clacluates the classificaiton output for each spotlight feature (sphere around a ith voxel, with predefined radius'''
    ''' iFeature is the index of Spotlight from AllPoints (embeded in DataArgDic)'''
    '''DataArcDic is a dict that contains all the information required for classification, including classification pipeline, conditions, location of the zstat file for each condition, spotligh points'''
    myPoints = DataArgDict['AllPoints'][iFeature] 
    # myPoints = SpherePoints(np.array([30,30,15]) , radius) # sample good point

    [X, Y, Gr]=MVPA_FeatureExtraction(DataArgDict['ZstatInfoSubj'],DataArgDict['CLconditions'],DataArgDict['Alldata'], myPoints) # extract data for classifer
    SpotPerformance=MVPA_Classification(X,Y,Gr, DataArgDict['pipeSpotLight'])
    return SpotPerformance 

    
def multi_run_wrapperLoopiFeatureDict(args):
    return LoopiFeatureDict(*args)    