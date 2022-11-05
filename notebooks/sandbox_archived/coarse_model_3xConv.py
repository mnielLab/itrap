#!/usr/bin/env python2.7

"""
NAME:        ProgramName
AUTHOR:      Helle Rus Povlsen 
DESCRIPTION:
	This is a template script.
"""

import sys
from sys import argv
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import h5py
import json

import tensorflow as tf
import keras

from keras.layers import Input
from keras.layers import Reshape
from keras.layers import Embedding
from keras.layers import Concatenate
from keras.layers import Conv2D
from keras.layers import MaxPooling1D
from keras.layers import MaxPooling2D
from keras.layers import AveragePooling2D
from keras.layers import Flatten
from keras.layers import BatchNormalization
from keras.layers import Dense
from keras.models import Model
from keras.callbacks import Callback
from keras.callbacks import ModelCheckpoint, EarlyStopping, TensorBoard
from keras.optimizers import SGD, RMSprop, Adagrad, Adadelta, Adam, Adamax, Nadam

# Improve GPU performance
gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction = 0.3, allow_growth=True, allocator_type='BFC')

from keras.backend.tensorflow_backend import set_session
set_session(tf.Session(config = tf.ConfigProto(gpu_options = gpu_options)))

# Custom activation function
from keras.layers import Activation
from keras import backend as K
from keras.utils.generic_utils import get_custom_objects

def swish_activation(z):
    return K.sigmoid(z) * z

get_custom_objects().update({'swish': Activation(swish_activation)})

# Evaluate model
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import accuracy_score
from scipy import interp

def ParseArguments():
    """ Handles program parameters and returns an argument class containing 
    all parameters """

    parser = ArgumentParser(description='Description of ProgramName')

    parser.add_argument('-g', '--grid_dir', metavar='DIR', help='Path to input file.',
            dest="grid_dir", required=True)
    parser.add_argument('-i', '--in_dir', metavar='DIR', help='Path to input files.',
            dest='in_dir', required=True)
    parser.add_argument('-o', '--out_dir', metavar='DIR', help='Path to output files.',
            dest='out_dir', required=True)
    parser.add_argument('-b', '--tb_dir', metavar='DIR', help='Path to output tensorboard files.',
            dest='tb_dir', required=True)
    parser.add_argument('-m', '--model_id', metavar='NAME', help='The run id in the grid search',
            dest='model_id', required=True)
    parser.add_argument('-s', '--setting', metavar='VAR', help='The setting for hyper-parameter search',
            dest='setting', nargs='?', default='coarse')
    parser.add_argument("-t", '--threads', metavar="N", type=int, dest="threads", 
            default=1, help="Number of processors to run in parallel [1].", required=False) 
    parser.add_argument('-v', '--verbose', action='store_true', dest="verbose", help="Verbose")
    
    args = parser.parse_args()

    return args


def coarse_parameters(grid_dir, model_name):
    import random
    from itertools import product

    bs = [128, 1024, 2048]
    lr = [0.001]
    of = ['RMSprop', 'Adagrad', 'Adam', 'Adam', 'Adam', 'Nadam']
    cn = [0, 1.0]
    af = ['relu', 'tanh', 'sigmoid', 'swish']
    po = ['max']
    bn = [True, False]
    do = [0.0]
    cf = [(10,10,10), (30,30,30)]
    nh = [16, 64]

    hyper_parameters = dict(batch_size = bs,
                            conv_activation_func = af,
                            conv_filters = cf,
                            norm = bn,

                            dense_dim = nh,
                            dense_activation_func = af,

                            optimizer = of,
                            clipnorm = cn)


    grid = list(product(*hyper_parameters.values()))

    rand_smpl = grid[random.sample(range(len(grid)), 1)[0]]

    sample_params = dict()
    for x,y in zip(hyper_parameters.keys(), rand_smpl):
        sample_params[x] = y

    with open(grid_dir + model_name + '.json', 'w') as fh:
        json.dump(sample_params, fh)

    return sample_params


def refine_parameters(grid_dir, model_name):
    import random
    from itertools import product

    bs = [1024] 
    lr = [0.001]
    of = ['Adam']
    cn = [0]
    ac = ['tanh']
    ad = ['swish']
    po = ['max','avg']
    bn = [False]
    do = [0.0]
    cf = [(30,30,30), (50,30,10), (30,20,10)]
    nh = [16, 32]
    mf = [[( 10,10), (10,1), (  1,2), ( 1,2), 'same'],
          [( 10, 5), (10,1), (  1,2), ( 1,2), 'same'],
          [( 10,10), ( 5,1), (  1,2), ( 1,2), 'same'],
          [( 10, 5), ( 5,1), (  1,2), ( 1,2), 'same'],
          [( 10,10), ( 2,1), (  5,2), ( 1,2), 'valid'],
          [( 10, 5), ( 2,1), (  5,2), ( 1,2), 'valid'],
          [(  5,10), ( 5,1), (  1,2), ( 1,2), 'same'],
          [(  5, 5), ( 5,1), (  1,2), ( 1,2), 'same'],
          [(  1,10), ( 1,1), (100,2), ( 1,2), 'same'],
          [(  1,10), ( 1,1), (100,2), ( 1,2), 'valid'],
          [(  1, 5), ( 1,1), (100,2), ( 1,2), 'same'],
          [(  1, 5), ( 1,1), (100,2), ( 1,2), 'valid'],
          [(  1,10), ( 1,1), ( 50,2), (50,2), 'same'],
          [(  1,10), ( 1,1), ( 50,2), (50,2), 'valid'],
          [(  1, 5), ( 1,1), ( 50,2), (50,2), 'same'],
          [(  1, 5), ( 1,1), ( 50,2), (50,2), 'valid'],
          [(  1,10), ( 1,1), ( 50,2), (25,2), 'same'],
          [(  1,10), ( 1,1), ( 50,2), (25,2), 'valid'],
          [(  1, 5), ( 1,1), ( 50,2), (25,2), 'same'],
          [(  1, 5), ( 1,1), ( 50,2), (25,2), 'valid'],
          [(100,10), ( 1,1), (  1,2), ( 1,2), 'valid'],
          [(100, 5), ( 1,1), (  1,2), ( 1,2), 'valid'],
          [(  1,21), ( 1,1), (  2,1), ( 2,1), 'valid'],
          [( 10,21), ( 1,1), (  2,1), ( 2,1), 'valid']]
 

    hyper_parameters = dict(batch_size = bs,
                            conv_activation_func = ac,
                            conv_pooling = po,
                            conv_filters = cf,
                            filters = mf,

                            dense_dim = nh,
                            dense_activation_func = ad)


    grid = list(product(*hyper_parameters.values()))

    rand_smpl = grid[random.sample(range(len(grid)), 1)[0]]

    sample_params = dict()
    for x,y in zip(hyper_parameters.keys(), rand_smpl):
        sample_params[x] = y

    print(sample_params)
    with open(grid_dir + model_name + '.json', 'w') as fh:
        json.dump(sample_params, fh)


    return sample_params


def generate_partitions(N_batches, n_folds=3):
    from sklearn.model_selection import KFold
    
    cv_indices = {'train':list(), 'test':list()}
    skf = KFold(n_splits=n_folds, shuffle=False, random_state=0)
    for train_index, test_index in skf.split(np.arange(N_batches)):
        cv_indices['train'].append(list(train_index))
        cv_indices['test'].append(list(test_index))
    
    return cv_indices


def load_batch(batch_index, input_directory):
    import h5py
    import random
    
    f = h5py.File(input_directory + str(batch_index) + ".hdf5", "r")
       
    images = f['images'][:]
    target = f['target'][:]
    depths = f['depths'][:]
        
    f.close()
    
    return ([images[f] for f in range(len(images))], np.array([[int(gt) for gt in genotype] for genotype in target]), depths)



def batch_iter(input_directory, CV_batches):
    import random
    
    while True:
        random.shuffle(CV_batches)
        
        for batch_index in CV_batches:
            yield load_batch(batch_index, input_directory)


def load_testset(test_indices, input_directory, batch_size):
    
    n_obsrv = len(test_indices)*batch_size
    n_reads = 100
    n_bases = 21
    n_genot = 10
    
    X_base = np.empty((n_obsrv, n_reads, n_bases), dtype=np.float32)
    X_baql = np.empty((n_obsrv, n_reads, n_bases), dtype=np.float32)
    X_mapq = np.empty((n_obsrv, n_reads, n_bases), dtype=np.float32)
    X_strd = np.empty((n_obsrv, n_reads, n_bases), dtype=np.float32)
    
    X_test = [X_base, X_baql, X_mapq, X_strd]
    y_test = np.empty((n_obsrv, n_genot), dtype=np.uint32)
    d_test = np.empty((n_obsrv), dtype=np.uint32)
    
    i, j = 0, batch_size

    for batch in test_indices:
        Xt, yt, dt = load_batch(batch, input_directory)
        
        X_test[0][i:j] = Xt[0]
        X_test[1][i:j] = Xt[1]
        X_test[2][i:j] = Xt[2]
        X_test[3][i:j] = Xt[3]
        y_test[i:j]    = yt
        d_test[i:j]    = dt
        
        i += batch_size
        j += batch_size
        
    return X_test, y_test, d_test



def get_batch_count(input_directory):
    import os
    
    return len([1 for f in os.listdir(input_directory) if f.endswith('.hdf5')]) - 2


def load_hyper_parameters(json_filename):
    with open(json_filename, "r") as fh:
        return json.load(fh)


def compile_model(shape_base = (100,21),
                  shape_qual = (100,21),
                  shape_mapq = (100,21),
                  shape_strd = (100,21),
                  shape_embd = (100,21),
                  
                  conv_filters = [30,30,30],
                  filters = [( 1, 5), ( 1,1), (  2,2), ( 2,2), 'valid'],
                  
                  conv_activation_func = 'tanh',
                  conv_pooling = 'max',
                  
                  dense_dim = 32,
                  dense_activation_func = 'swish',

                  drop = 0,
                  norm = False,
                  
                  output_dim = 10,
                  output_activation_func = 'softmax',
                  
                  loss = 'categorical_crossentropy',
                  optimizer = 'adam',
                  learning_rate = 0.001,
                  clipnorm = 0,
                  metrics = ['accuracy']):

    #assert type(input_shape)==tuple,"Input shape must be a tuple"
    
    ######################################################################################################
    ############################################## Initiate ##############################################
    ######################################################################################################

    if filters[0][0] == 100:
        conv_filter_dim1 = filters[0]
        conv_filter_dim2 = (1,5)
    elif filters[0][1] == 21:
        conv_filter_dim1 = filters[0]
        conv_filter_dim2 = (5,1)
    else:
        conv_filter_dim1 = filters[0]
        conv_filter_dim2 = filters[0]

    if filters[0][1] == 10 and filters[4] == 'valid':
        conv_filter_dim1 = filters[0]
        conv_filter_dim2 = (1,5)

    conv_strides = filters[1]
    conv_padding = filters[-1]
    
    conv_drop = drop
    conv_norm = norm

    pool_dim = filters[2] 
    pool_stride = filters[3]

    dense_drop = drop
    dense_norm = norm

    ######################################################################################################
    ############################################ Define input ############################################
    ######################################################################################################
    l_base = Input(shape=shape_base, name='input_base', dtype=np.uint8)
    l_qual = Input(shape=shape_qual, name='input_qual')
    l_mapq = Input(shape=shape_mapq, name='input_mapq')
    l_strd = Input(shape=shape_strd, name='input_strd', dtype=np.uint8)

    ######################################################################################################
    ########################################### Reshape input ############################################
    ######################################################################################################
    l_qual_reshape = Reshape(target_shape=shape_qual + (1,), name='qual_reshape')(l_qual)
    l_mapq_reshape = Reshape(target_shape=shape_mapq + (1,), name='mapq_reshape')(l_mapq)

    ######################################################################################################
    ############################################ Embed input #############################################
    ######################################################################################################
    l_embed_base = Embedding(input_dim=6, output_dim=4, mask_zero=False, input_length=shape_embd, name='base_embed')(l_base)
    l_embed_strd = Embedding(input_dim=3, output_dim=2, mask_zero=False, input_length=shape_embd, name='strd_embed')(l_strd)

    ######################################################################################################
    ######################################## Concatenate features ########################################
    ######################################################################################################
    l_main = Concatenate(axis=3)([l_embed_base, l_qual_reshape, l_mapq_reshape, l_embed_strd])

    ######################################################################################################
    ######################################### Convolutional layer ########################################
    ######################################################################################################
    l_conv = Conv2D(filters=conv_filters[0], 
                    kernel_size=conv_filter_dim1, 
                    strides=conv_strides, 
                    padding=conv_padding, 
                    activation=conv_activation_func)(l_main)
    
    l_conv = Conv2D(filters=conv_filters[1], 
                    kernel_size=conv_filter_dim2, 
                    strides=conv_strides, 
                    padding=conv_padding, 
                    activation=conv_activation_func,)(l_conv)
    
    l_conv = Conv2D(filters=conv_filters[2], 
                    kernel_size=conv_filter_dim2, 
                    strides=conv_strides, 
                    padding=conv_padding, 
                    activation=conv_activation_func,)(l_conv)
    
    ######################################################################################################
    ############################################    Pooling    ###########################################
    ######################################################################################################
    if conv_pooling == 'max':
        l_pool = MaxPooling2D(pool_size=pool_dim, strides=pool_stride)(l_conv)
    else:
        l_pool = AveragePooling2D(pool_size=pool_dim, strides=pool_stride)(l_conv)
        
    if conv_drop:
        l_pool = Dropout(conv_drop)(l_pool)
    if conv_norm:
        l_pool = BatchNormalization()(l_pool)
        
    l_flat = Flatten()(l_pool)
    
    ######################################################################################################
    ############################################## Dense layer ###########################################
    ######################################################################################################
    l_dens = Dense(units=dense_dim, activation=dense_activation_func)(l_flat)
    
    if dense_drop:
        l_dens = Dropout(dense_drop)(l_dens)
    if dense_norm:
        l_dens = BatchNormalization()(l_dens)

    ######################################################################################################
    ############################################# Output layer ###########################################
    ######################################################################################################
    l_out = Dense(units=output_dim, activation=output_activation_func, name='output')(l_dens)

    ######################################################################################################
    ############################################ Compile model ###########################################
    ######################################################################################################
    model = Model(inputs=[l_base, l_qual, l_mapq, l_strd], outputs=[l_out])
    
    optimizer_opts = dict()
    if learning_rate:
        optimizer_opts.update(lr=learning_rate)
    #if decay:
    #    optimizer_opts.update(decay=decay)
    #if clipvalue:
    #    optimizer_opts.update(clipvalue=clipvalue)
    if clipnorm:
        optimizer_opts.update(clipnorm=clipnorm)
    
    if optimizer == 'SGD':
        optimizer = SGD(**optimizer_opts) # momentum?
    elif optimizer == 'RMSprop':
        optimizer = RMSprop(**optimizer_opts)
    elif optimizer == 'Adagrad':
        optimizer = Adagrad(**optimizer_opts)
    elif optimizer == 'Adadelta':
        optimizer = Adadelta(**optimizer_opts)
    elif optimizer == 'Adam':
        optimizer = Adam(**optimizer_opts)
    elif optimizer == 'Adamax':
        optimizer = Adamax(**optimizer_opts)
    elif optimizer == 'Nadam':
        optimizer = Nadam(**optimizer_opts)
    
    
    model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
    
    #model.predict_proba = model.predict  what is this?!
    
    return model


class IntervalEvaluation(Callback):
    """
    Evaluate model at each epoch end using implemented metrics. 
    Stop training when a monitored quantity has stopped improving.
    """

    def __init__(self, X, y, model_name, interval=1, patience=20):
        super(Callback, self).__init__()

        self.interval = interval
        self.X_test, self.y_test = X,y
        self.auc = []
        self.pr = []
        self.mcc = []
        self.model_name = model_name #"/home/hrp/Master/data/input_test/model_weights/%s.h5"
        self.burnin = 4 #30
        self.optimal_epoch = 0
        
        self.monitor = ['auc', 'pr', 'mcc']
        self.monitor_value = None ### OBS! self.auc, self.pr, self.mcc ?!
        self.patience = patience
        self.wait = 0
        self.stopped_epoch = 0
        
    def get_metrics(self, epoch):
        y_pred = self.model.predict(self.X_test)
        
        if 'auc' in self.monitor:
            self.auc.append( roc_auc_score(self.y_test, y_pred, average="micro") )
            print("AUC %.2f" %self.auc[epoch])
        if 'pr' in self.monitor:
            self.pr.append( average_precision_score( self.y_test, y_pred, average="micro" ) )
            print("PR %.2f" %self.pr[epoch])
        if 'mcc' in self.monitor:
            y_pred_cm = y_pred.argmax(1)
            y_test_cm = self.y_test.argmax(1)
            self.mcc.append( matthews_corrcoef(y_test_cm, y_pred_cm) )
            print("MCC %.2f" %self.mcc[epoch])
            
        try:
            self.monitor_value = np.array(self.auc) + np.array(self.pr) + np.array(self.mcc)
        except ValueError as error:
            print("Check which monitor metrics are used. ", error)

        
    def on_train_begin(self, logs=None): # When is this useful?
        # Allow instances to be re-used
        self.wait = 0
        self.stopped_epoch = 0

    def on_epoch_end(self, epoch, logs={}):
        if epoch % self.interval == 0:
            self.get_metrics(epoch)
            if epoch > self.burnin and self.auc[epoch] > 0.5:
                if np.argmax( self.monitor_value ) == epoch: # This method finds plateaus -> hence, stops before decline
                    self.wait = 0
                    self.optimal_epoch = epoch
                    self.model.save_weights(self.model_name)
                else:
                    self.wait += 1
                    if self.wait >= self.patience:
                        self.stopped_epoch = epoch
                        self.model.stop_training = True

    def on_train_end(self, logs={}):
        self.model.load_weights(self.model_name)
        
        fh = h5py.File(self.model_name, "a")
        performance = fh.create_group("performance")
        performance.create_dataset("auc", data=self.auc)
        performance.create_dataset("pr", data=self.pr)
        performance.create_dataset("mcc", data=self.mcc)
        performance.attrs['stopped_epoch'] = self.stopped_epoch + 1
        performance.attrs['optimal_epoch'] = self.optimal_epoch + 1
        performance.attrs['optimal_auc'] = self.auc[self.optimal_epoch]
        performance.attrs['optimal_pr'] = self.pr[self.optimal_epoch]
        performance.attrs['optimal_mcc'] = self.mcc[self.optimal_epoch]
        fh.close()
        
        if self.stopped_epoch > 0:
            print('Epoch %d: optimal performance' % (self.optimal_epoch + 1))
            print('Epoch %d: early stopping' % (self.stopped_epoch + 1))
    
    


def main(argv):	
    args = ParseArguments()

    ###################
    #   Initialize    #
    ###################

    if args.setting == 'coarse':
        sample_hyperparams = coarse_parameters(args.grid_dir, args.model_id)
    else:
        sample_hyperparams = refine_parameters(args.grid_dir, args.model_id)

    batch_size = sample_hyperparams.pop('batch_size', None)

    input_dir = args.in_dir + '/%s/' % batch_size
    
    N_batches = get_batch_count(input_dir)

    cv_indices = generate_partitions(N_batches)


    #####################
    #   Compile model   #
    #####################
    try:
        model = compile_model(**sample_hyperparams) #**sample_hyperparams
    except Exception as error:
        print("Network couldn't compile with the hyper-parameters given in id: %s" %args.model_id)
        print(error)
        sys.exit(1)

    trainable_parameters = model.count_params()
    if trainable_parameters > 500000:
        print("Too many parameters: %s" %trainable_parameters)
        print(model.summary())
        sys.exit(1)
    else:
        print(trainable_parameters)

    ####################
    #   Train model    #
    ####################
    for fold in range(len(cv_indices['train'])): #len(partitions['train'])
        model = compile_model(**sample_hyperparams) #**sample_hyperparams
        model_name = "%s/%s_%s.h5" %(args.out_dir, args.model_id, fold)
        
        X_test, y_test, d_test = load_testset(cv_indices['test'][fold],
                                              input_directory=input_dir,
                                              batch_size=batch_size)

        
        #alt_performance_metrics = IntervalEvaluation(X_test, y_test, model_name, interval=1)

        early_stopping = EarlyStopping(monitor='val_acc',
                                       mode='max',
                                       patience=20, 
                                       min_delta=0.001,
                                       verbose=1)

        model_checkpoint = ModelCheckpoint(model_name,
                                           monitor='val_acc',
                                           mode='max',
                                           verbose=0,
                                           save_best_only=True,
                                           save_weights_only=True)

        tensorboard = TensorBoard(log_dir=args.tb_dir + args.model_id + "_" + str(fold),  
                                  write_images=True,
                                  histogram_freq=0)


        # Fit the model
        history = model.fit_generator(batch_iter(input_directory=input_dir, 
                                                 CV_batches=cv_indices['train'][fold]), 
                                      steps_per_epoch=60, #len(cv_indices['train'][fold]), 
                                      epochs=200, 
                                      callbacks=[early_stopping, model_checkpoint, tensorboard],
                                      validation_data=(X_test, y_test),
                                      verbose=2)
        
        # Predict
        pred = model.predict(X_test, verbose=0)
        
        fh = h5py.File(model_name, "a")
        performance = fh.create_group("performance")
        performance.create_dataset("acc", data=history.history['acc'])
        performance.create_dataset("loss", data=history.history['loss'])
        performance.create_dataset("val_acc", data=history.history['val_acc'])
        performance.create_dataset("val_loss", data=history.history['val_loss'])
        performance.attrs['stopped_epoch'] = early_stopping.stopped_epoch + 1
        performance.attrs['optimal_epoch'] = early_stopping.stopped_epoch - early_stopping.patience + 1
        performance.attrs['optimal_acc'] = history.history['val_acc'][early_stopping.stopped_epoch - early_stopping.patience]
        performance.attrs['optimal_loss'] = history.history['val_loss'][early_stopping.stopped_epoch - early_stopping.patience]

        prediction = fh.create_group("prediction")
        prediction.create_dataset("prediction", data=pred)
        prediction.create_dataset("y_test", data=y_test)
        fh.close()
        

    print("Successful training completed")
	
if __name__ == '__main__':
   main(argv)
