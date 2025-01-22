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
os.chdir("Y:/data/Feature_selection_framework/final_test_all_wilc/")

def df2array(df, YellowPage):
    valid_uid = YellowPage.loc[YellowPage.iloc[:,2] == 'Test']['uid']
    train_uid = YellowPage.loc[YellowPage.iloc[:,2] != 'Test']['uid']
    valid_sub = df.loc[valid_uid,]
    train_sub = df.loc[train_uid,]
    #make sure all data is float
    valid_values = valid_sub.values
    train_values = train_sub.values
    #split the X and Y in train and valid
    X = (train_values[:,:train_values.shape[1]-1].astype(np.float32), valid_values[:,:valid_values.shape[1]-1].astype(np.float32))# X of train and test
    Y = (train_values[:,train_values.shape[1]-1].astype(np.intc), valid_values[:,valid_values.shape[1]-1].astype(np.intc))#Y of train and test
    return X, Y

def MLP_model(Data, Yellow_Page, parm): #parm = best hp
    #split the train
    x, y = df2array(Data, Yellow_Page) # use first fold for reference

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

def MAT_PROCESS (MAT_NAME, VAR, INPUT, CV_OUTPUT, YELLOW_PAGE):
    #get the best hp
    cv_out_sub = CV_OUTPUT[CV_OUTPUT['mat'] == MAT_NAME] #subset by the mat name
    perform_list = [] #place holder
    for key in VAR:
        #key = 'filter'
        data_pv = INPUT[VAR[key]]
        best_parm = cv_out_sub.loc[cv_out_sub['Feature'] == key, ["n_neuron", "learning_rate"]].values.flatten().tolist()

        data_pv['genotype'] = [i.split("_")[2] for i in data_pv.index] #parse the genotype from the uid
        #encode the genotype
        encoder = LabelEncoder()
        data_pv['genotype'] = encoder.fit_transform(data_pv['genotype']) #0 is Mut, 1 is WT

        #build the MLP model
        test_df = MLP_model(data_pv, YELLOW_PAGE, best_parm) #dont need index
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

    #read the cv output for the mat
    cv_out_file = '../CV_wilc/MLP/MLP_' + f_name + '.csv' #get CV output
    cv_out = pandas.read_csv(cv_out_file, sep = ',', decimal = '.')
    return(MAT_PROCESS(mat_name, var_list, input_data, cv_out, yellow_page))

f_name = "Opt_on"
#get the corresponding yellow page filename
yellow_page_prefix = '../int_output/'
yellow_page_suffix = '_kFolds_yellow_page.csv'
yellow_page_file = yellow_page_prefix + f_name + yellow_page_suffix
#read the yellow page
yellow_page = pandas.read_csv(yellow_page_file, sep = ",", decimal='.')

on_files = glob.glob("../transformed_data/*/*on.csv")

try:
   mp.set_start_method('spawn', force=True)
   print("spawned")
except RuntimeError:
   pass

if __name__ == '__main__':
    with mp.Pool(5) as p:
        main_list = p.map(main, on_files)
    MAT_auc_df = pandas.concat(main_list) #concatenate dfs
    MAT_auc_df.insert(loc = 0, column = 'light', value = 'on') #add the light condition
    MAT_auc_df.insert(loc = 0, column = 'classifer', value = 'MLP') #add the classifier

    #set the output directory and filenames
    out_path = "./MLP/"
    if not os.path.exists(out_path):
        os.makedirs(out_path)#create the output dir for each dataset
    #save the accuracy dataframe to an output csv file
    MAT_auc_df.to_csv(out_path + "MLP_" + f_name + '_test_wilc.csv', sep=',', index = False)
