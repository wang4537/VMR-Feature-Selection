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
os.chdir("Y:/data/Feature_selection_framework/6HV_preprocess_sep/")

def df2array(Set1, Set2):
    #make sure all data is float
    Set1_values = Set1.values
    Set2_values = Set2.values
    #split the X and Y in train and valid
    X = (Set1_values[:,:Set1_values.shape[1]-1].astype(np.float32), Set2_values[:,:Set2_values.shape[1]-1].astype(np.float32))# X of train and test
    Y = (Set1_values[:,Set1_values.shape[1]-1].astype(np.intc), Set2_values[:,Set2_values.shape[1]-1].astype(np.intc))#Y of train and test
    return X, Y

def MLP_model(TRAIN, TEST, parm): #parm = best hp
    #split the train
    x, y = df2array(TRAIN, TEST) # use first fold for reference

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

def MAT_PROCESS (MAT_NAME, VAR, INPUT_TRAIN, INPUT_TEST, CV_OUTPUT):
    #get the best hp
    cv_out_sub = CV_OUTPUT[CV_OUTPUT['mat'] == MAT_NAME] #subset by the mat name
    perform_list = [] #place holder
    for key in VAR:
        #key = 'filter'
        train_pv = INPUT_TRAIN[VAR[key]]
        test_pv = INPUT_TEST[VAR[key]]
        best_parm = cv_out_sub.loc[cv_out_sub['Feature'] == key, ["n_neuron", "learning_rate"]].values.flatten().tolist()

        train_pv['genotype'] = [i.split("_")[2] for i in train_pv.index] #parse the genotype from the uid
        test_pv['genotype'] = [i.split("_")[2] for i in test_pv.index] #parse the genotype from the uid
        #encode the genotype
        encoder = LabelEncoder()
        train_pv['genotype'] = encoder.fit_transform(train_pv['genotype']) #0 is Mut, 1 is WT
        test_pv['genotype'] = encoder.fit_transform(test_pv['genotype']) #0 is Mut, 1 is WT

        #build the MLP model
        test_df = MLP_model(train_pv, test_pv, best_parm) #dont need index
        #add column names
        best_parm_df = pandas.DataFrame(best_parm, index=["n_neuron", "learning_rate"]).T
        out_df = pandas.concat([best_parm_df, test_df], axis=1)
        #add column for the feature name
        out_df.insert(loc = 0, column = 'Feature', value = key)
        perform_list.append(out_df)

    perform_df = pandas.concat(perform_list) #concatenatate dfs
    perform_df.insert(loc = 0, column = 'mat', value = MAT_NAME) #add matricization name
    print('MAT complete:' + MAT_NAME)
    return(perform_df)

def main(x):
    prefix = '../3transformed_data_sep/'#head dir name
    mat_name = x[len(prefix):].split('\\')[0] #remove dir and file extension
    input_train = pandas.read_csv(x, sep = ",", decimal='.', index_col=0)
    test_file = prefix + mat_name + '/' + mat_name + '_' + f_name + '_test.csv'
    input_test =pandas.read_csv(test_file, sep = ",", decimal='.', index_col=0)

    filter_file = glob.glob('../4feature_selection_sep/filter_wilc/' + mat_name + '/*off_vol.csv', recursive=True) #get filter features
    filter_list = pandas.read_csv(filter_file[0])['feature']

    embedded_file = glob.glob('../4feature_selection_sep/embedded_wilc/' + mat_name + '/*off_sigVar_pass_list.csv', recursive=True) #get the embedded features
    embedded_list = pandas.read_csv(embedded_file[0]).iloc[:,0]
    intersect_list = list(set(filter_list) & set(embedded_list)) #intersect
    union_list = sorted(list(set(filter_list)) + list(set(embedded_list))) #union

    #list of feature sets
    var_list = {'full': input_train.columns[1:]}
    if(len(filter_list) != 0):
        var_list['filter'] = filter_list
    if(len(embedded_list) != 0):
        var_list['embedded'] = embedded_list
    if(len(intersect_list) != 0):
        var_list['intersect'] = intersect_list
    if(len(union_list) !=0):
        var_list['union'] = union_list

    #read the cv output for the mat
    cv_out_file = '../5CV_sep/MLP/MLP_' + f_name + '_CV.csv' #get CV output
    cv_out = pandas.read_csv(cv_out_file, sep = ',', decimal = '.')
    return(MAT_PROCESS(mat_name, var_list, input_train, input_test, cv_out))

f_name = "Opt_off"
off_files = glob.glob("../3transformed_data_sep/*/*off_train.csv")

try:
   mp.set_start_method('spawn', force=True)
   print("spawned")
except RuntimeError:
   pass

if __name__ == '__main__':
    with mp.Pool(3) as p:
        main_list = p.map(main, off_files)
    MAT_auc_df = pandas.concat(main_list) #concatenate dfs
    MAT_auc_df.insert(loc = 0, column = 'light', value = 'off') #add the light condition
    MAT_auc_df.insert(loc = 0, column = 'classifer', value = 'MLP') #add the classifier

    #set the output directory and filenames
    out_path = "./MLP/"
    if not os.path.exists(out_path):
        os.makedirs(out_path)#create the output dir for each dataset
    #save the accuracy dataframe to an output csv file
    MAT_auc_df.to_csv(out_path + "MLP_" + f_name + '_test_wilc.csv', sep=',', index = False)
