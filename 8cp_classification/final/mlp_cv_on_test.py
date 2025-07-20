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
os.chdir("/depot/yleung/data/Feature_selection_framework/8cp_classification/final/")

def MLP_model(x, y, parm): #parm = best hp
    #design the model
    model = Sequential([
        keras.layers.Dense(int(parm[0]), input_shape = (x[0].shape[1],), activation = 'relu'),
        keras.layers.Dense(int(parm[0]), activation = 'relu'),
        keras.layers.Dense(1, activation = 'sigmoid')])
    model.compile(loss = keras.losses.BinaryCrossentropy(from_logits=True),
                  optimizer = keras.optimizers.Adam(learning_rate=parm[1]),
                  metrics = ['accuracy'])
    #fit the model
    history = model.fit(x[0], y[0], epochs = 200, batch_size = 64, verbose = 0) #training X and training Y

    #predict the validating set
    valid_pred = model.predict(x[1])
    #confusion matrix
    valid_matrix = metrics.confusion_matrix(y[1], np.rint(valid_pred))
    #break down
    tn = valid_matrix[0, 0]
    tp = valid_matrix[1, 1]
    fn = valid_matrix[1, 0]
    fp = valid_matrix[0, 1]

    #output metrics
    auroc = metrics.roc_auc_score(y[1], np.rint(valid_pred)) #area under ROC curve
    acc = metrics.accuracy_score(y[1], np.rint(valid_pred)) #accuracy
    sns = metrics.recall_score(y[1], np.rint(valid_pred)) #sensitivity = recall
    sps = tn / (tn + fp) #specificity
    prc = metrics.precision_score(y[1], np.rint(valid_pred)) #precision
    kappa = metrics.cohen_kappa_score(y[1], np.rint(valid_pred)) #Cohen's kappa

    perform_metrics = pandas.DataFrame([auroc, acc, sns, sps, prc, kappa], index= ['AUROC', 'Accuracy', 'Sensitivity', 'Specificity','Precision', 'Kappa']).T
    print('Test performance calculated!')
    return  perform_metrics

# MAIN
f_name = "Opt_on"

#training data
train =  pandas.read_csv("../project_data/Opt_on_cp_train_data.csv")
#testing data
test = pandas.read_csv("../project_data/Opt_on_cp_test_data.csv")
#encode the genotype
encoder = LabelEncoder()
train['genotype'] = encoder.fit_transform(train['genotype']) #0 is Mut, 1 is WT
test['genotype'] = encoder.fit_transform(test['genotype']) #0 is Mut, 1 is WT
#read the best parameter
cv_out = pandas.read_csv("../CV/MLP/MLP_Opt_on_CP_cv.csv", sep = ',', decimal = '.')
best_parm = cv_out[["n_neuron", "learning_rate"]].values.flatten().tolist()
#make sure all data is float
test_values = test.values
train_values = train.values
#split the X and Y in train and valid
X = (train_values[:,2:].astype(np.float32), test_values[:,2:].astype(np.float32))# X of train and test
Y = (train_values[:,1].astype(np.intc), test_values[:,1].astype(np.intc))#Y of train and test
#run MLP
out_df = MLP_model(X, Y, best_parm)
#add columns to output
out_df = pandas.concat([cv_out[["n_neuron", "learning_rate"]], out_df], axis=1)
out_df.insert(loc = 0, column = 'Feature', value = "default")
out_df.insert(loc = 0, column = 'mat', value = "CP") #add matricization name
out_df.insert(loc = 0, column = 'light', value = 'on') #add the light condition
out_df.insert(loc = 0, column = 'classifer', value = 'MLP') #add the classifier
#save output
out_path = "./MLP/"
if not os.path.exists(out_path):
    os.makedirs(out_path)#create the output dir for each dataset
out_df.to_csv(out_path + "MLP_" + f_name + '_CP_test.csv', sep=',', index = False)
