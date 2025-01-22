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
os.chdir("Y:/data/Feature_selection_framework/CV_wilc/")

def df2array(df, YellowPage, fold):
    valid_uid = YellowPage.loc[YellowPage.iloc[:,fold+1] == 'Valid']['uid']
    train_uid = YellowPage.loc[YellowPage.iloc[:,fold+1] == 'Train']['uid']
    valid_sub = df.loc[valid_uid,]
    train_sub = df.loc[train_uid,]
    #make sure all data is float
    valid_values = valid_sub.values
    train_values = train_sub.values
    #split the X and Y in train and valid
    X = (train_values[:,:train_values.shape[1]-1].astype(np.float32), valid_values[:,:valid_values.shape[1]-1].astype(np.float32))# X of train and test
    Y = (train_values[:,train_values.shape[1]-1].astype(np.intc), valid_values[:,valid_values.shape[1]-1].astype(np.intc))#Y of train and test
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
    history = model.fit(x[0], y[0], epochs = 200, batch_size = 64, verbose = 0) #training X and training Y

    #predict the validating set
    valid_pred = model.predict(x[1])
    #confusion matrix
    valid_matrix = metrics.confusion_matrix(y[1], np.rint(valid_pred))
    #output metrics
    auc = metrics.roc_auc_score(y[1], np.rint(valid_pred)) #area under ROC curve
 #AUROC
    return  auc

#make a function to load, preprocess and build the MLP model
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

def MAT_PROCESS (MAT_NAME, VAR, INPUT):
    max_auc_list = []
    for key in VAR:
        parm_acc_list = []
        data_pv = INPUT[VAR[key]] #subset the data by column names
        data_pv['genotype'] = [i.split("_")[2] for i in data_pv.index] #parse the genotype from the uid
        #encode the genotype
        encoder = LabelEncoder()
        data_pv['genotype'] = encoder.fit_transform(data_pv['genotype']) #0 is Mut, 1 is WT
        #build the MLP model
        parm_acc_list.append(CV(data_pv, yellow_page, MLP_parm_list))
        #convert the list to data frame
        MLP_parm_df = pandas.DataFrame(MLP_parm_list, columns=["n_neuron", "learning_rate"])
        #add the
        out_df = pandas.DataFrame(parm_acc_list, index=['AUROC'])
        #add column names
        out_df = pandas.concat([MLP_parm_df, out_df.T], axis=1)
        #find the max value
        out_df_mx = out_df[out_df.AUROC == out_df.AUROC.max()]
        #add column for the feature name
        out_df_mx.insert(loc = 0, column = 'Feature', value = key)
        max_auc_list.append(out_df_mx)

    max_auc_df = pandas.concat(max_auc_list) #concatenatate dfs
    max_auc_df.insert(loc = 0, column = 'mat', value = MAT_NAME) #add matricization name
    print('MAT complete:' + MAT_NAME)
    return(max_auc_df)

def main(x):
    prefix = '../transformed_data/'#head dir name
    mat_name = x[len(prefix):].split('\\')[0] #remove dir and file extension
    input_data = pandas.read_csv(x, sep = ",", decimal='.', index_col=0)

    filter_file = glob.glob('../filter_wilc/' + mat_name + '/*on_vol.csv', recursive=True) #get filter features
    filter_list = pandas.read_csv(filter_file[0])['feature']

    embedded_file = glob.glob('../embedded_wilc/' + mat_name + '/*on_sigVar_pass_list.csv', recursive=True) #get the embedded features
    embedded_list = pandas.read_csv(embedded_file[0]).iloc[:,0]
    intersect_list = list(set(filter_list) & set(embedded_list)) #intersect
    union_list = sorted(list(set(filter_list)) + list(set(embedded_list))) #union

    #list of feature sets
    var_list = {'full': input_data.columns[1:]}
    if(len(filter_list) != 0):
        var_list['filter'] = filter_list
    if(len(embedded_list) != 0):
        var_list['embedded'] = embedded_list
    if(len(intersect_list) != 0):
        var_list['intersect'] = intersect_list
    if(len(union_list) !=0):
        var_list['union'] = union_list

    return(MAT_PROCESS(mat_name, var_list, input_data))

f_name = "Opt_on"
#get the corresponding yellow page filename
yellow_page_prefix = '../int_output/'
yellow_page_suffix = '_kFolds_yellow_page.csv'
yellow_page_file = yellow_page_prefix + f_name + yellow_page_suffix
#read the yellow page
yellow_page = pandas.read_csv(yellow_page_file, sep = ",", decimal='.')

#MLP hyperparameters
MLP_parm_list = [[x,y] for x in [8,16,32,64,128] for y in [0.0001, 0.001, 0.01]] #x = number of neurons, y = learning rate for 'adam'

on_files = glob.glob("../transformed_data/*/*on.csv")

try:
   mp.set_start_method('spawn', force=True)
   print("spawned")
except RuntimeError:
   pass


if __name__ == '__main__':
    with mp.Pool(7) as p:
        main_list = p.map(main, on_files)

    MAT_auc_df = pandas.concat(main_list) #concatenate dfs
    MAT_auc_df.insert(loc = 0, column = 'light', value = 'on') #add the light condition
    MAT_auc_df.insert(loc = 0, column = 'classifer', value = 'MLP') #add the classifier

    #set the output directory and filenames
    out_path = "./MLP/"
    if not os.path.exists(out_path):
        os.makedirs(out_path)#create the output dir for each dataset
    #save the accuracy dataframe to an output csv file
    MAT_auc_df.to_csv(out_path + "MLP_" + f_name + '.csv', sep=',', index = False)
