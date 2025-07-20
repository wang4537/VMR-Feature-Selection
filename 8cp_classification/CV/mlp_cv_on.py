#load the packages
#import numpy
import pandas
#import io
import os
import glob
from sklearn.preprocessing import LabelEncoder
import numpy as np
from sklearn import metrics
#following required for LSTM
from tensorflow import keras
from keras.models import Sequential
#from keras.layers import LSTM
#from keras.layers import Dropout
#import keras_tuner as kt
import tensorflow as tf
import multiprocessing as mp

tf.get_logger().setLevel('ERROR')
#set the working directory
os.chdir("/depot/yleung/data/Feature_selection_framework/8cp_classification/CV/")

def df2array(df, YellowPage, fold):
    valid_uid = YellowPage.loc[YellowPage.iloc[:,fold+1] == 'Valid']['uid']
    train_uid = YellowPage.loc[YellowPage.iloc[:,fold+1] == 'Train']['uid']
    valid_sub = df.loc[valid_uid,]
    train_sub = df.loc[train_uid,]
    #make sure all data is float
    valid_values = valid_sub.values
    train_values = train_sub.values
    #split the X and Y in train and valid
    X = (train_values[:,1:].astype(np.float32), valid_values[:,1:].astype(np.float32))# X of train and test
    Y = (train_values[:,0].astype(np.intc), valid_values[:,0].astype(np.intc))#Y of train and test
    return X, Y

def MLP_model(x, y, parm):
    #design the model
    model = Sequential([
        keras.layers.Dense(parm[0], input_shape = (x[0].shape[1],), activation = 'relu'),
        keras.layers.Dense(parm[0], activation = 'relu'),
        keras.layers.Dense(1, activation = 'sigmoid')])
    model.compile(loss = keras.losses.BinaryCrossentropy(from_logits=True),
                  optimizer = keras.optimizers.Adam(learning_rate=parm[1]),
                  metrics = ['accuracy'])
    #fit the model
    history = model.fit(x[0], y[0], epochs = 200, batch_size = 8, verbose = 0) #training X and training Y

    #predict the validating set
    valid_pred = model.predict(x[1])
    #confusion matrix
    valid_matrix = metrics.confusion_matrix(y[1], np.rint(valid_pred))
    print(valid_matrix)
    #output metrics
    auc = metrics.roc_auc_score(y[1], np.rint(valid_pred)) #area under ROC curve
 #AUROC
    return  auc

def CV(Data, Yellow_Page, MLP_parm):
    #create a list to store the accuracy for each fold
    MLP_fold_perm = []
    for k in range(1,11):
        #split the train
        X, Y = df2array(Data, Yellow_Page, k)
        #apply each model for the independent and response variables
        k_acc =[]
        for i in MLP_parm:
            k_acc.append(MLP_model(X, Y, i))
        MLP_fold_perm.append(k_acc)
        #print an indicator
        print("Fold Complete:", k)
    #calculate the mean and folds
    cv_mean = np.mean(MLP_fold_perm, axis = 0)
    print('Average CV calculated!')
    return(cv_mean)
#### MAIN ####
f_name = "Opt_on"
#get the corresponding yellow page filename
yellow_page_prefix = '../../2preprocess_sep/'
yellow_page_suffix = '_kFolds_yellow_page.csv'
yellow_page_file = yellow_page_prefix + f_name + yellow_page_suffix
#read the yellow page
yellow_page = pandas.read_csv(yellow_page_file, sep = ",", decimal='.')

#MLP hyperparameters
MLP_parm_list = [[x,y] for x in [8,16,32,64,128] for y in [0.0001, 0.001, 0.01]] #x = number of neurons, y = learning rate for 'adam'
input_file = "../project_data/Opt_on_cp_train_data.csv"
input_data = pandas.read_csv(input_file, sep = ",", decimal='.', index_col=0)
#encode the genotype
encoder = LabelEncoder()
input_data['genotype'] = encoder.fit_transform(input_data['genotype']) #0 is Mut, 1 is WT
#build the MLP model
parm_acc_list = CV(input_data, yellow_page, MLP_parm_list)

out_df = pandas.DataFrame(parm_acc_list, columns=['AUROC'])
#convert the list to data frame
MLP_parm_df = pandas.DataFrame(MLP_parm_list, columns=["n_neuron", "learning_rate"])
#add column names
out_df = pandas.concat([MLP_parm_df, out_df], axis=1)
#find the max value
out_df_mx = out_df[out_df.AUROC == out_df.AUROC.max()]
#add column for the feature name
out_df_mx.insert(loc = 0, column = 'Feature', value = "default")
out_df_mx.insert(loc = 0, column = 'mat', value = "CP") #add matricization name
out_df_mx.insert(loc = 0, column = 'light', value = 'on') #add the light condition
out_df_mx.insert(loc = 0, column = 'classifer', value = 'MLP') #add the classifier
#set the output directory and filenames
out_path = "./MLP/"
if not os.path.exists(out_path):
    os.makedirs(out_path)#create the output dir for each dataset
#save the accuracy dataframe to an output csv file
out_df_mx.to_csv(out_path + "MLP_" + f_name + '_CP_cv.csv', sep=',', index = False)
