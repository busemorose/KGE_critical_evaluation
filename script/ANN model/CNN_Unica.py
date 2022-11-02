# -*- coding: utf-8 -*-
"""
@author: Andreas Wunsch
andreas.wunsch@kit.edu
wunsch.andreas.edu@gmail.com

"""
#reproducability
from numpy.random import seed
seed(1)
import tensorflow as tf
tf.random.set_seed(1)

import numpy as np
from bayes_opt import BayesianOptimization
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs
import os
import pandas as pd
import datetime
from scipy import stats,optimize
from matplotlib import pyplot
from sklearn.preprocessing import MinMaxScaler
from keras.utils.vis_utils import plot_model
import json
import tensorflow as tf


gpus = tf.config.experimental.list_physical_devices('GPU')

#%% functions
def load_data():
    filepath = "./unica_ANN_dataset_dummy.txt"
    data = pd.read_csv(filepath,parse_dates=['date'],index_col=0, decimal = '.', sep=';')
    
    # print(data.isna().sum())
    # data['S_postojna'] = data['S_postojna'] .interpolate(method='linear',limit = 1)
    # data['rH'] = data['rH'] .interpolate(method='linear',limit = 1)
    # data['PET'] = data['PET'] .interpolate(method='linear')
    # print(data.index[(data['PET'].isna())])
    
    #add sinus signal
    tt = np.linspace(0, data.shape[0]-1, data.shape[0])
    T = np.asarray(data['T'])
    res = fit_sin(tt, T)
    data['Tsin'] = res["fitfunc"](tt)
    
    #drop unnecessary data
    data.drop(columns=['P_thiessen','SR_precip'],inplace=True)
    
    return data

def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = np.array(tt)
    yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}


def split_data(data, GLOBAL_SETTINGS):
    #separate data and create overlaps for best data usage of input data (_ext)
    dataset = data[(data.index < GLOBAL_SETTINGS["test_start"])] 
    dataset1 = data[(data.index < GLOBAL_SETTINGS["optset_start"])]
    
    TrainingData = dataset[(dataset.index < GLOBAL_SETTINGS["stopset_start"])]
    
    StopData = dataset[(dataset.index >= GLOBAL_SETTINGS["stopset_start"]) & (dataset.index < GLOBAL_SETTINGS["optset_start"])]
    StopData_ext = pd.concat([TrainingData.iloc[-GLOBAL_SETTINGS["seq_length"]:], StopData], axis=0)
    
    OptData = dataset[(dataset.index >= GLOBAL_SETTINGS["optset_start"]) & (dataset.index < GLOBAL_SETTINGS["test_start"])]
    OptData_ext = pd.concat([dataset1.iloc[-GLOBAL_SETTINGS["seq_length"]:], OptData], axis=0)

    TestData = data[(data.index >= GLOBAL_SETTINGS["test_start"]) & (data.index <= GLOBAL_SETTINGS["test_end"])] 
    TestData_ext = pd.concat([dataset.iloc[-GLOBAL_SETTINGS["seq_length"]:], TestData], axis=0) 
    
    return TrainingData, StopData, StopData_ext, OptData, OptData_ext, TestData, TestData_ext

def to_supervised(data, GLOBAL_SETTINGS):
    # reformat data in sequence format
    X, Y = list(), list()
    # step over the entire history one time step at a time
    for i in range(len(data)):
        # find the end of this pattern
        end_idx = i + GLOBAL_SETTINGS["seq_length"]
        # check if we are beyond the dataset
        if end_idx >= len(data):
            break
        # gather input and output parts of the pattern
        seq_x, seq_y = data[i:end_idx, 1:], data[end_idx, 0]
        X.append(seq_x)
        Y.append(seq_y)
    return np.array(X), np.array(Y)

class MCDropout(tf.keras.layers.Dropout):
    #define Monte Carlo Dropout Layer
    def call(self, inputs):
        return super().call(inputs, training=True)
    
def Qmodel(ini,GLOBAL_SETTINGS,X_train):
    # define model
    seed(ini+37657)
    tf.random.set_seed(ini+37657)
    
    inp = tf.keras.Input(shape=(GLOBAL_SETTINGS["seq_length"], X_train.shape[2]))
    cnn = tf.keras.layers.Conv1D(filters=GLOBAL_SETTINGS["filters"],
                                         kernel_size=GLOBAL_SETTINGS["kernel_size"],
                                         activation='relu')(inp)
    cnn = tf.keras.layers.MaxPool1D(padding='same')(cnn)
    cnn = MCDropout(GLOBAL_SETTINGS["dropout"])(cnn)
    
    cnn = tf.keras.layers.Flatten()(cnn)
    cnn = tf.keras.layers.Dense(GLOBAL_SETTINGS["dense_size"], activation='relu')(cnn)
    output1 = tf.keras.layers.Dense(1, activation='linear')(cnn)

    # tie together
    model = tf.keras.Model(inputs=inp, outputs=output1)
    optimizer = tf.keras.optimizers.Adam(learning_rate=GLOBAL_SETTINGS["learning_rate"],
                                         epsilon=10E-3, clipnorm=GLOBAL_SETTINGS["clip_norm"])
    
    model.compile(loss='mse',optimizer=optimizer, metrics=['mse'])
    
    return model

def predict_distribution(X, model, n):
    # preds = [model(X, training=True) for _ in range(n)]
    preds = [model(X) for _ in range(n)]
    return np.hstack(preds)

def simulate_testset(densesize, seqlength, batchsize, filters,PET,T,rH,S_postojna,nS_postojna,S_cerknica,nS_cerknica,Tsin):
    
    GLOBAL_SETTINGS['batch_size'] = batchsize
    GLOBAL_SETTINGS['dense_size'] = densesize
    GLOBAL_SETTINGS['filters'] = filters
    GLOBAL_SETTINGS['seq_length'] = seqlength
    
    ## load data
    data = load_data()
    
    if not PET:
        data.drop(columns=['PET'],inplace=True)
    if not rH:
        data.drop(columns=['rH'],inplace=True)
    if not T:
        data.drop(columns=['T'],inplace=True)
    if not S_postojna:
        data.drop(columns=['S_postojna'],inplace=True)
    if not nS_postojna:
        data.drop(columns=['nS_postojna'],inplace=True)
    if not S_cerknica:
        data.drop(columns=['S_cerknica'],inplace=True)
    if not nS_cerknica:
        data.drop(columns=['nS_cerknica'],inplace=True)
    if not Tsin:
        data.drop(columns=['Tsin'],inplace=True)
    
    GS = json.dumps(str(GLOBAL_SETTINGS))
    with open("./global_settings.json", "w") as outfile: 
        outfile.write(GS) 
    
        
    #split and scale data
    
    TrainingData, StopData, StopData_ext, OptData, OptData_ext, TestData, TestData_ext = split_data(data, GLOBAL_SETTINGS)
        
    scaler_Q = MinMaxScaler(feature_range=(0, 1))
    scaler_Q.fit(pd.DataFrame(TrainingData[['Qobs']]))
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(pd.DataFrame(TrainingData))
    data_n = pd.DataFrame(scaler.transform(data), index=data.index, columns=data.columns)

    TrainingData_n, StopData_n, StopData_ext_n, OptData__n, OptData_ext_n, TestData_n, TestData_ext_n = split_data(data_n, GLOBAL_SETTINGS)
    
    #sequence data
    X_train, Y_train = to_supervised(TrainingData_n.values, GLOBAL_SETTINGS)
    X_stop, Y_stop = to_supervised(StopData_ext_n.values, GLOBAL_SETTINGS)
    X_opt, Y_opt = to_supervised(StopData_ext_n.values, GLOBAL_SETTINGS) 
    X_test, Y_test = to_supervised(TestData_ext_n.values, GLOBAL_SETTINGS) 

    #build and train model with different initializations
    inimax = 10
    n_distr = 100 #number of predictions per ini (MCDropout)

    testresults_members_all = np.zeros((len(X_test), inimax*n_distr))

    #loop initializations
    for ini in range(inimax):
        model = Qmodel(ini,GLOBAL_SETTINGS,X_train)
        
        # early stopping
        es = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=0, 
                                              patience=GLOBAL_SETTINGS["patience"],restore_best_weights = True)
    
        # fit network
        history = model.fit(x=X_train,y=Y_train, validation_data=(X_stop,Y_stop), 
                            epochs=GLOBAL_SETTINGS["epochs"], verbose=GLOBAL_SETTINGS["verbose"],
                            batch_size=GLOBAL_SETTINGS["batch_size"], callbacks=[es])

        # plot loss during training
        pyplot.figure(figsize=(10,4))
        pyplot.title('Loss')
        pyplot.plot(history.history['loss'], label='train')
        pyplot.plot(history.history['val_loss'], label='val_loss')
        pyplot.ylabel('Loss', size=12)
        pyplot.xlabel('Epochs',size=12)
        pyplot.legend()
        pyplot.savefig('Train_history_ini_'+str(ini)+'.png', dpi=300)
        pyplot.show()
        
        # save model
        model_name = 'model_ini' + str(ini)
        model.save('./' + model_name)
        
        y_pred_distribution = predict_distribution(X_test, model, n_distr)
        test_sim = scaler_Q.inverse_transform(y_pred_distribution)
        testresults_members_all[:,ini*n_distr:ini*n_distr+n_distr]=test_sim
    
        pyplot.plot(scaler_Q.inverse_transform(Y_test.reshape(-1,1)),'k')
        pyplot.plot(test_sim.mean(axis=1),'r',alpha = 0.5)
        pyplot.show()

    sim1_uncertainty = [np.quantile(testresults_members_all, 0.05, axis=1),np.quantile(testresults_members_all, 0.95, axis=1)]
    
    plot_model(model, to_file='model_plot.png', show_shapes=True, show_layer_names=True, dpi=300)
    
    test_sim_mean1 = np.mean(testresults_members_all,axis = 1)    
    sim1 = np.asarray(test_sim_mean1.reshape(-1,1))
    
    
    Y_test_n = Y_test
    Y_test = scaler_Q.inverse_transform(Y_test_n.reshape(-1,1))
    obs1 = Y_test.reshape(-1,1)

    err = sim1-obs1
    err_rel = (sim1-obs1)/(np.max(data['Qobs'])-np.min(data['Qobs']))
    err_nash = obs1 - np.mean(obs1)

    NSE = 1 - ((np.sum(err ** 2)) / (np.sum((err_nash) ** 2)))
    try:
        r = stats.pearsonr(sim1[:,0], obs1[:,0])
        r = r[0] #r
    except:
        r = [np.nan, np.nan]
        r = r[0] #r
    
    R2 = r ** 2
    
    RMSE =  np.sqrt(np.mean(err ** 2))
    rRMSE = np.sqrt(np.mean(err_rel ** 2)) * 100
    Bias = np.mean(err)
    rBias = np.mean(err_rel) * 100
    # PI = np.nan #1 - ((np.sum(err ** 2)) / (np.sum((err_PI) ** 2)))
    
    alpha = np.std(sim1)/np.std(obs1)
    beta = np.mean(sim1)/np.mean(obs1)
    KGE = 1-np.sqrt((r-1)**2+(alpha-1)**2+(beta-1)**2) #KGE
    
    scores = pd.DataFrame(np.array([[NSE, R2, RMSE, rRMSE, Bias, rBias, KGE]]),
                   columns=['NSE','R2','RMSE','rRMSE','Bias','rBias','KGE'])
    print(scores)
    

    return scores, TestData, sim1, obs1, inimax, testresults_members_all, sim1_uncertainty

def bayesOpt_function(densesize, seqlength, batchsize, filters,PET,T,rH,S_postojna,nS_postojna,S_cerknica,nS_cerknica,Tsin):
    
    #translate into usabable parameters
    densesize = 2**int(densesize)
    seqlength = int(seqlength)
    batchsize = 2**int(batchsize)
    filters = 2**int(filters)
    
    PET = int(round(PET))
    T = int(round(T))
    rH = int(round(rH))
    S_postojna = int(round(S_postojna))
    nS_postojna = int(round(nS_postojna))
    S_cerknica = int(round(S_cerknica))
    nS_cerknica = int(round(nS_cerknica))
    Tsin = int(round(Tsin))
    
    #store in dict
    GLOBAL_SETTINGS['batch_size'] = batchsize
    GLOBAL_SETTINGS['dense_size'] = densesize
    GLOBAL_SETTINGS['filters'] = filters
    GLOBAL_SETTINGS['seq_length'] = seqlength
    
    ## load data
    data = load_data()
    
    #drop data according to optimization
    if not PET:
        data.drop(columns=['PET'],inplace=True)
    if not rH:
        data.drop(columns=['rH'],inplace=True)
    if not T:
        data.drop(columns=['T'],inplace=True)
    if not S_postojna:
        data.drop(columns=['S_postojna'],inplace=True)
    if not nS_postojna:
        data.drop(columns=['nS_postojna'],inplace=True)
    if not S_cerknica:
        data.drop(columns=['S_cerknica'],inplace=True)
    if not nS_cerknica:
        data.drop(columns=['nS_cerknica'],inplace=True)
    if not Tsin:
        data.drop(columns=['Tsin'],inplace=True)
        

    #split and scale data
    TrainingData, StopData, StopData_ext, OptData, OptData_ext, TestData, TestData_ext = split_data(data, GLOBAL_SETTINGS)
        
    scaler_Q = MinMaxScaler(feature_range=(0, 1))
    scaler_Q.fit(pd.DataFrame(TrainingData[['Qobs']]))
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaler.fit(pd.DataFrame(TrainingData))
    data_n = pd.DataFrame(scaler.transform(data), index=data.index, columns=data.columns)

    TrainingData_n, StopData_n, StopData_ext_n, OptData__n, OptData_ext_n, TestData_n, TestData_ext_n = split_data(data_n, GLOBAL_SETTINGS)
    
    #sequence data
    X_train, Y_train = to_supervised(TrainingData_n.values, GLOBAL_SETTINGS)
    X_stop, Y_stop = to_supervised(StopData_ext_n.values, GLOBAL_SETTINGS)
    X_opt, Y_opt = to_supervised(OptData_ext_n.values, GLOBAL_SETTINGS) 
    # X_test, Y_test = to_supervised(TestData_ext_n.values, GLOBAL_SETTINGS) 

    #build and train model with different initializations
    inimax = 1
    optresults_members = np.zeros((len(X_opt), inimax))

    for ini in range(inimax):
        
        model = Qmodel(ini,GLOBAL_SETTINGS,X_train)
        
        # early stopping
        es = tf.keras.callbacks.EarlyStopping(monitor='val_loss', mode='min', verbose=0, 
                                              patience=GLOBAL_SETTINGS["patience"],restore_best_weights = True)
    
        # fit network
        history = model.fit(x=X_train,y=Y_train, validation_data=(X_stop,Y_stop), 
                            epochs=GLOBAL_SETTINGS["epochs"], verbose=GLOBAL_SETTINGS["verbose"],
                            batch_size=GLOBAL_SETTINGS["batch_size"], callbacks=[es]) 
        
        # plot loss during training
        pyplot.figure(figsize=(10,4))
        pyplot.title('Loss')
        pyplot.plot(history.history['loss'], label='train')
        pyplot.plot(history.history['val_loss'], label='val_loss')
        pyplot.ylabel('Loss', size=12)
        pyplot.xlabel('Epochs',size=12)
        pyplot.legend()
        pyplot.show()
        
        opt_sim_n = model.predict(X_opt)
        opt_sim = scaler_Q.inverse_transform(opt_sim_n)
        optresults_members[:, ini] = opt_sim[:,0].reshape(-1,)
        
        
    opt_sim_mean1 = np.mean(optresults_members,axis = 1)    
    sim1 = np.asarray(opt_sim_mean1.reshape(-1,1))

    Y_opt_n = Y_opt
    Y_opt = scaler_Q.inverse_transform(Y_opt_n.reshape(-1,1))
    obs1 = Y_opt.reshape(-1,1)
    
    pyplot.plot(obs1,'k')
    pyplot.plot(sim1,'b',alpha = 0.5)
    pyplot.show()

    # get scores
    err = sim1-obs1
    MSE = np.mean(err ** 2)

    return (-1)*MSE

class newJSONLogger(JSONLogger) :

      def __init__(self, path):
            self._path=None
            super(JSONLogger, self).__init__()
            self._path = path if path[-5:] == ".json" else path + ".json"
            


#%% start optimization

with tf.device("/gpu:0"): # or cpu
    
    time1 = datetime.datetime.now()
    basedir = './'
    os.chdir(basedir)


    time_single = datetime.datetime.now()
    seed(1+37657)
    tf.random.set_seed(1+37657)
    
    # Bounded region of parameter space
    pbounds = {'seqlength': (2,365), 
               'densesize': (4,9),
               'batchsize': (4,9),
               'filters': (4,9),
               'PET': (0,1),
               'T': (0,1),
               'rH': (0,1),
               'S_postojna': (0,1),
               'nS_postojna': (0,1),
               'S_cerknica': (0,1),
               'nS_cerknica': (0,1),
               'Tsin': (0,1)}

    #define counters for optimization
    optsteps1 = 30 # random initial steps
    optsteps2 = 60 # least no of steps
    optsteps3 = 15 # how many steps no improvement
    optsteps4 = 150 # max no of steps
    
    GLOBAL_SETTINGS = {
        'kernel_size': 3, 
        'dropout': 0.1,
        'clip_norm': True,
        'epochs': 150,
        'patience': 20,
        'learning_rate': 1e-3,
        'verbose': 1,
        'stopset_start': pd.to_datetime('01102006', format='%d%m%Y'),
        'optset_start': pd.to_datetime('01102008', format='%d%m%Y'),
        'test_start': pd.to_datetime('01102010', format='%d%m%Y'),
        'test_end': pd.to_datetime('31122018', format='%d%m%Y')
    }
    
    #HP optimizer
    optimizer = BayesianOptimization(
        f= bayesOpt_function, #fct to optimize
        pbounds=pbounds, #HP boundaries to optimize within
        random_state=1, 
        verbose = 0 # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent, verbose = 2 prints everything
        )
       
    # #load existing optimizer
    log_already_available = 0
    if os.path.isfile("./logs.json"):
        load_logs(optimizer, logs=["./logs.json"]);
        print("\nExisting optimizer is already aware of {} points.".format(len(optimizer.space)))
        log_already_available = 1
    
    # Save progress
    logger = newJSONLogger(path="./logs.json")
    optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)
    
    if log_already_available == 0:
        optimizer.maximize(
                init_points=optsteps1, #steps of random exploration 
                n_iter=0, # steps of bayesian optimization
                acq="ei",# ei  = expected improvmenet (probably the most common acquisition function) 
                xi=0.05  #  Prefer exploitation (xi=0.0) / Prefer exploration (xi=0.1)
                )
    
    # optimize while improvement 
    current_step = len(optimizer.res)
    beststep = False
    step = -1
    while not beststep:
        step = step + 1
        beststep = optimizer.res[step] == optimizer.max #aktuell beste Iteration suchen

    while current_step < optsteps2: 
            current_step = len(optimizer.res)
            beststep = False
            step = -1
            while not beststep:
                step = step + 1
                beststep = optimizer.res[step] == optimizer.max
            print("\nbeststep {}, current step {}".format(step+1, current_step+1))
            optimizer.maximize(
                init_points=0, #steps of random exploration 
                n_iter=1, # steps of bayesian optimization
                acq="ei",# ei  = expected improvmenet (probably the most common acquisition function) 
                xi=0.05  #  Prefer exploitation (xi=0.0) / Prefer exploration (xi=0.1)
                )
            
    while (step + optsteps3 > current_step and current_step < optsteps4): 
            current_step = len(optimizer.res)
            beststep = False
            step = -1
            while not beststep:
                step = step + 1
                beststep = optimizer.res[step] == optimizer.max
                
            print("\nbeststep {}, current step {}".format(step+1, current_step+1))
            optimizer.maximize(
                init_points=0, #steps of random exploration 
                n_iter=1, # steps of bayesian optimization
                acq="ei",# ei  = expected improvmenet (probably the most common acquisition function) 
                xi=0.05  #  Prefer exploitation (xi=0.0) / Prefer exploration (xi=0.1)
                )
        
    print("\nBEST:\t{}".format(optimizer.max))
        
    #get best values from optimizer
    densesize = 2**int(optimizer.max.get("params").get("densesize"))
    seqlength = int(optimizer.max.get("params").get("seqlength"))
    batchsize = 2**int(optimizer.max.get("params").get("batchsize"))
    filters = 2**int(optimizer.max.get("params").get("filters"))
    
    PET = int(round(optimizer.max.get("params").get("PET")))
    T = int(round(optimizer.max.get("params").get("T")))
    rH = int(round(optimizer.max.get("params").get("rH")))
    S_postojna = int(round(optimizer.max.get("params").get("S_postojna")))
    nS_postojna = int(round(optimizer.max.get("params").get("nS_postojna")))
    S_cerknica = int(round(optimizer.max.get("params").get("S_cerknica")))
    nS_cerknica = int(round(optimizer.max.get("params").get("nS_cerknica")))
    Tsin = int(round(optimizer.max.get("params").get("Tsin")))


    #run test set simulations

    scores, TestData, sim1, obs1, inimax, testresults_members_all, sim1_uncertainty = simulate_testset(densesize, seqlength, batchsize, filters,PET,T,rH,S_postojna,nS_postojna,S_cerknica,nS_cerknica,Tsin)

#%% plot

lb = sim1_uncertainty[0]
ub = sim1_uncertainty[1]
sim = sim1
obs = obs1

pyplot.figure(figsize=(10,3))
pyplot.fill_between(TestData.index, lb,
                    ub, facecolor = (1,0.7,0,0.99),
                    label ='90% confidence',linewidth = 0.8,
                    edgecolor = (1,0.7,0,0.99))

pyplot.plot(TestData.index, sim, color = 'r', label ="simulated mean", alpha=0.8,linewidth=0.8)
pyplot.plot(TestData.index, obs, 'k', label ="observed", linewidth=0.7,alpha=0.5)
pyplot.title("Unica Springs", size=15)
pyplot.ylabel('Q [m³/s]', size=12)
pyplot.xlabel('Date',size=12)
pyplot.legend(fancybox = False, framealpha = 0, edgecolor = 'k',loc='upper right')
pyplot.grid(b=True, which='major', color='#666666', alpha = 0.1, linestyle='-')
pyplot.tight_layout()

s = """NSE\nR²\nRMSE\nBias\nKGE"""
pyplot.figtext(0.08, 0.6, s)

s = """{:.2f}\n{:.2f}\n{:.2f}\n{:.2f}\n{:.2f}\n""".format(scores.NSE[0],scores.R2[0],
scores.RMSE[0],scores.Bias[0],scores.KGE[0])
pyplot.figtext(0.13, 0.55, s)
    
pyplot.savefig('Plot_Unica.png', dpi=500)
pyplot.show()


#%% save results
printdf = pd.DataFrame(data=np.c_[obs,sim,lb,ub],index=TestData.index)
printdf = printdf.rename(columns={0: 'Obs', 1: 'Sim', 2:'lb:0.05', 3:'ub:95'})
printdf.to_csv('./results.txt',sep=';', float_format = '%.6f')