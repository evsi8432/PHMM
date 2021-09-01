import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import sys

from copy import deepcopy

from scipy.stats import gamma
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde
from scipy.special import expit
from scipy.special import logit
from scipy.special import logsumexp
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.signal import convolve
from scipy.interpolate import interp1d

import divebomb

import time
import pickle

import importlib

import Preprocessor
import Parameters
import HHMM
import Visualisor

rand_seed = int(sys.argv[1])
model = str(sys.argv[2])

random.seed(rand_seed)
np.random.seed(rand_seed)


### initialize parameters ####

ndives = 100
ndatasets = 1

# dive duration parameters
dd_mu = np.array([20.0,80.0])
dd_sig = np.array([5.0,50.0])

dd_shape = np.square(dd_mu)/np.square(dd_sig)
dd_scale = np.square(dd_sig)/np.array(dd_mu)

# FoVeDBA parameters (this is for the SQUARE of the FoVeDBA)
FoVeDBA_sin_shape = np.array([[(5.0/n) for n in range(1,52)],[(5.0/n) for n in range(1,52)]])
FoVeDBA_sin_scale = np.ones((2,51))
FoVeDBA_sin_shape[1,1] = 250
FoVeDBA_sin_shape[1,2] = 50

# average acceleration parameters
acc_mu = np.array([0.0,0.0])
acc_sig = np.array([0.05,0.1])

# number of states for each substate
K0 = 2
K1 = 2

# randomly initialize a probablity transition matrix
ptm_crude = np.array([[0.5,0.5],
                      [0.5,0.5]])

ptm_fine = [np.array([[0.5,0.5],
                      [0.1,0.9]]),
            np.array([[0.8,0.2],
                      [0.3,0.7]])]

# randomly intialize a correlation within states
corr_crude = [0.0,0.0]
corr_fine = [0.99,0.95]

# train new models?
train_new = True

# initialize the initial states
delta_crude = np.ones(K0)/K0
for _ in range(100):
    delta_crude = delta_crude.dot(ptm_crude)

delta_fines = []
for k0 in range(K0):
    delta_fine = np.ones(K1)/K1
    for _ in range(100):
        delta_fine = delta_fine.dot(ptm_fine[k0])
    delta_fines.append(delta_fine)


### Create Data ###

def create_data():

    data = []
    data_V = []
    data_FV = []
    freqs = np.fft.rfftfreq(100, d=1/50.0)
    thresh = 5
    thresh_ind = max(np.where(freqs <= thresh)[0]) + 1

    time = 0

    for dive_num in range(ndives):

        # select dive type
        if dive_num == 0:
            dive_type = np.random.choice(K0,p=delta_crude)
            dd_mu_t = np.copy(dd_mu[dive_type])
        else:
            dive_type_tm1 = dive_type
            dive_type = np.random.choice(K0,p=ptm_crude[dive_type,:])
            dd_mu_t = np.copy(dd_mu[dive_type])
            dd_mu_t = corr_crude[dive_type]*dd_tm1 + (1.0-corr_crude[dive_type])*dd_mu_t

        # select dive duration
        dd_sig_t = dd_sig[dive_type]
        dd_shape_t = np.square(dd_mu_t)/np.square(dd_sig_t)
        dd_scale_t = np.square(dd_sig_t)/np.array(dd_mu_t)

        dive_duration = gamma.rvs(dd_shape_t,0,dd_scale_t)
        dd_tm1 = dive_duration
        datum = {'dive_type': dive_type, 'dive_duration': dive_duration}
        datum_V = {'dive_type': dive_type, 'dive_duration': dive_duration}
        datum_FV = {'dive_type': dive_type, 'dive_duration': dive_duration}
        nsegs = int(dive_duration/2.0)

        subdive_features = []
        subdive_features_V = []
        subdive_features_FV = []

        for seg_num in range(nsegs):

            seg = {}
            seg_V = [{},{}]
            seg_FV = {}

            # find seg type
            if seg_num == 0:
                subdive_type = np.random.choice(K1,p=delta_fines[dive_type])
                FoVeDBA_sin_mu_t = FoVeDBA_sin_shape[subdive_type,:]*FoVeDBA_sin_scale[subdive_type,:]
                acc_mu_t = np.copy(acc_mu[subdive_type])
            else:
                subdive_type = np.random.choice(K1,p=ptm_fine[dive_type][subdive_type,:])
                FoVeDBA_sin_mu_t = FoVeDBA_sin_shape[subdive_type,:]*FoVeDBA_sin_scale[subdive_type,:]
                acc_mu_t = np.copy(acc_mu[subdive_type])
                acc_mu_t = corr_fine[subdive_type]*acc_tm1 + (1.0-corr_fine[subdive_type])*acc_mu_t

            seg['subdive_type'] = subdive_type
            seg_FV['subdive_type'] = subdive_type
            seg_V[0]['subdive_type'] = subdive_type
            seg_V[1]['subdive_type'] = subdive_type

            # find average acceleration
            acc_sig_t = acc_sig[subdive_type]

            # find FoVeDBA
            FoVeDBA_sin_sig2_t = FoVeDBA_sin_shape[subdive_type,:]*FoVeDBA_sin_scale[subdive_type,:]**2
            FoVeDBA_sin_shape_t = np.square(FoVeDBA_sin_mu_t)/FoVeDBA_sin_sig2_t
            FoVeDBA_sin_scale_t = FoVeDBA_sin_sig2_t/FoVeDBA_sin_mu_t
            FoVeDBA_sin = gamma.rvs(FoVeDBA_sin_shape_t,0,FoVeDBA_sin_scale_t)
            FoVeDBA_sin_tm1 = np.copy(FoVeDBA_sin)

            seg['FoVeDBA_full'] = FoVeDBA_sin.T
            seg['FoVeDBA'] = np.sum(FoVeDBA_sin.T[1:thresh_ind]) + 0.001
            seg_FV['FoVeDBA'] = np.sum(FoVeDBA_sin.T[1:thresh_ind]) + 0.001

            # find VeDBA
            pm = (2*np.random.choice(2,size = FoVeDBA_sin.shape)-1)
            FoVeDBA = pm*np.sqrt(FoVeDBA_sin)*1.0j
            A = np.array(np.fft.irfft(FoVeDBA))
            V0 = norm.rvs(acc_mu_t,acc_sig_t)
            acc_tm1 = V0
            A += np.linspace(V0,V0,100).T

            seg['A'] = A
            seg['avg_A'] = np.mean(A)
            seg_V[0]['A'] = np.mean(A[:50])
            seg_V[1]['A'] = np.mean(A[50:])
            seg_FV['A'] = np.mean(A)

            # find time
            seg['time'] = np.arange(time,time+2,1/50)
            time += 2

            subdive_features.append(seg)
            subdive_features_V.append(seg_V[0])
            subdive_features_V.append(seg_V[1])
            subdive_features_FV.append(seg_FV)

        datum['subdive_features'] = subdive_features
        datum_V['subdive_features'] = subdive_features_V
        datum_FV['subdive_features'] = subdive_features_FV

        data.append(datum)
        data_V.append(datum_V)
        data_FV.append(datum_FV)

    return data,data_V,data_FV


### Set initial Parameters ###

### CarHMM ###
hmm_FV_theta = [{'dive_duration': {'mu': np.array([np.mean(dd_mu)]),
                                   'sig': np.array([np.sqrt(np.mean(dd_sig**2) + np.var(dd_mu))]),
                                   'corr': np.array([-10.])}},
                 [{'FoVeDBA': {'mu': np.sum(FoVeDBA_sin_shape,1),
                               'sig': np.sqrt(np.sum(FoVeDBA_sin_shape,1)),
                               'corr': np.array([ -10., -10.])},
                   'A': {'mu': acc_mu,
                         'sig': acc_sig,
                         'corr': logit(corr_fine)}}]]

ptm_fine_temp = ptm_fine[0]/2 + ptm_fine[1]/2

hmm_FV_eta = [np.array([[0]]),
             [np.array([[0,np.log(1.0/ptm_fine_temp[0,0]-1)],
                        [np.log(1.0/ptm_fine_temp[1,1]-1),0]])]]


### CarHHMM, no Z2 ###
hhmm_V_theta = [{'dive_duration': {'mu': dd_mu,
                                   'sig': dd_sig,
                                   'corr': np.array([-10.,-10.])}},
                 [{'A': {'mu': acc_mu,
                         'sig': acc_sig,
                         'corr': logit(corr_fine)}},
                  {'A': {'mu': acc_mu,
                         'sig': acc_sig,
                         'corr': logit(corr_fine)}}]]

hhmm_V_eta = [np.array([[0,np.log(1.0/ptm_crude[0,0]-1)],
                        [np.log(1.0/ptm_crude[1,1]-1),0]]),

             [np.array([[0,np.log(1.0/ptm_fine[0][0,0]-1)],
                        [np.log(1.0/ptm_fine[0][1,1]-1),0]]),

              np.array([[0,np.log(1.0/ptm_fine[1][0,0]-1)],
                        [np.log(1.0/ptm_fine[1][1,1]-1),0]])]]


### HHMM ###
hhmm_FV_uncorr_theta = [{'dive_duration': {'mu': dd_mu,
                                           'sig': dd_sig,
                                           'corr': np.array([-10.,-10.])}},
                         [{'FoVeDBA': {'mu': np.sum(FoVeDBA_sin_shape,1),
                                       'sig': np.sqrt(np.sum(FoVeDBA_sin_shape,1)),
                                       'corr': np.array([ -10., -10.])},
                           'A': {'mu': acc_mu,
                                 'sig': acc_sig,
                                 'corr': np.array([-10.,-10.])}},
                          {'FoVeDBA': {'mu': np.sum(FoVeDBA_sin_shape,1),
                                       'sig': np.sqrt(np.sum(FoVeDBA_sin_shape,1)),
                                       'corr': np.array([ -10., -10.])},
                           'A': {'mu': acc_mu,
                                 'sig': acc_sig,
                                 'corr': np.array([-10.,-10.])}}]]

hhmm_FV_uncorr_eta = [np.array([[0,np.log(1.0/ptm_crude[0,0]-1)],
                                [np.log(1.0/ptm_crude[1,1]-1),0]]),

                      [np.array([[0,np.log(1.0/ptm_fine[0][0,0]-1)],
                                 [np.log(1.0/ptm_fine[0][1,1]-1),0]]),

                       np.array([[0,np.log(1.0/ptm_fine[1][0,0]-1)],
                                 [np.log(1.0/ptm_fine[1][1,1]-1),0]])]]

### CarHHMM ###
hhmm_FV_theta = [{'dive_duration': {'mu': dd_mu,
                                    'sig': dd_sig,
                                    'corr': np.array([-10.,-10.])}},
                 [{'FoVeDBA': {'mu': np.sum(FoVeDBA_sin_shape,1),
                               'sig': np.sqrt(np.sum(FoVeDBA_sin_shape,1)),
                               'corr': np.array([ -10., -10.])},
                   'A': {'mu': acc_mu,
                         'sig': acc_sig,
                         'corr': logit(corr_fine)}},
                 {'FoVeDBA': {'mu': np.sum(FoVeDBA_sin_shape,1),
                               'sig': np.sqrt(np.sum(FoVeDBA_sin_shape,1)),
                               'corr': np.array([ -10., -10.])},
                   'A': {'mu': acc_mu,
                         'sig': acc_sig,
                         'corr': logit(corr_fine)}}]]

hhmm_FV_eta = [np.array([[0,np.log(1.0/ptm_crude[0,0]-1)],
                         [np.log(1.0/ptm_crude[1,1]-1),0]]),

               [np.array([[0,np.log(1.0/ptm_fine[0][0,0]-1)],
                          [np.log(1.0/ptm_fine[0][1,1]-1),0]]),

                np.array([[0,np.log(1.0/ptm_fine[1][0,0]-1)],
                          [np.log(1.0/ptm_fine[1][1,1]-1),0]])]]


### Train Model ###

datasets = []
datasets_V = []
datasets_FV = []

hmm_FVs = []
hhmm_Vs = []
hhmm_FV_uncorrs = []
hhmm_FVs = []

h = 0.001

for dataset_num in range(ndatasets):


    ### generate data ###
    print('')
    print('GENERATING DATA')
    print('')
    data,data_V,data_FV = create_data()


    ### CarHMM ###
    if model == 'CarHMM':
        print('')
        print('STARTING CarHMM')
        print('')
        pars = Parameters.Parameters()
        pars.K = [1,2]
        pars.features = [{'dive_duration':{'corr':False,'f':'gamma'}},
                         {'FoVeDBA':{'corr':False,'f':'gamma'},
                          'A':{'corr':True,'f':'normal'}}]
        pars.theta = hmm_FV_theta

        hmm_FV = HHMM.HHMM(pars,data_FV)
        hmm_FV.theta = hmm_FV_theta
        hmm_FV.eta = hmm_FV_eta
        hmm_FV.true_theta = deepcopy(hmm_FV_theta)
        hmm_FV.true_eta = deepcopy(hmm_FV_eta)

        hmm_FV.train_DM(data_FV)
        hmm_FV.get_SEs(data_FV,h)
        data,data_V,data_FV = create_data()

        # make test data
        for dive_num,datum in enumerate(data_FV):
            _,_,posts,_ = hmm_FV.fwd_bwd(datum['subdive_features'],[1,0])
            for i,post in enumerate(posts.T):
                data[dive_num]['subdive_features'][i]['hmm_FV_dive'] = 0.0
                data[dive_num]['subdive_features'][i]['hmm_FV_subdive'] = post[1]

        # get confusion matrix
        CM_coarse = np.zeros((2,2))
        CM_fine = [np.zeros((2,2)),np.zeros((2,2))]
        for dive in data:
            dive_type = dive['dive_type']
            CM_coarse[dive_type,1] += dive['subdive_features'][0]['hmm_FV_dive']
            CM_coarse[dive_type,0] += max(0,1.0-dive['subdive_features'][0]['hmm_FV_dive'])
            for seg in dive['subdive_features']:
                subdive_type = seg['subdive_type']
                CM_fine[dive_type][subdive_type,1] += seg['hmm_FV_subdive']
                CM_fine[dive_type][subdive_type,0] += max(0,1.0-seg['hmm_FV_subdive'])

        hmm_FV.CM = [CM_coarse,CM_fine]

        # save data
        hmm_FV.save('../Params/hmm_FV_%d_%d'%(dataset_num,rand_seed))


    ### HHMM ###
    elif model == 'HHMM':
        print('')
        print('STARTING HHMM')
        print('')
        pars = Parameters.Parameters()
        pars.K = [2,2]
        pars.features = [{'dive_duration':{'corr':False,'f':'gamma'}},
                         {'FoVeDBA':{'corr':False,'f':'gamma'},
                          'A':{'corr':False,'f':'normal'}}]

        hhmm_FV_uncorr = HHMM.HHMM(pars,data_FV)
        hhmm_FV_uncorr.theta = hhmm_FV_uncorr_theta
        hhmm_FV_uncorr.eta = hhmm_FV_uncorr_eta
        hhmm_FV_uncorr.true_theta = deepcopy(hhmm_FV_uncorr_theta)
        hhmm_FV_uncorr.true_eta = deepcopy(hhmm_FV_uncorr_eta)

        hhmm_FV_uncorr.train_DM(data_FV)
        hhmm_FV_uncorr.get_SEs(data_FV,h)
        data,data_V,data_FV = create_data()

        # crude posterior
        _,_,posts_crude,_ = hhmm_FV_uncorr.fwd_bwd(data_FV,[0])

        # fine posterior
        for dive_num,datum in enumerate(data_FV):
            data[dive_num]['hhmm_FV_uncorr'] = posts_crude.T[dive_num,1]
            _,_,posts_fine_0,_ = hhmm_FV_uncorr.fwd_bwd(datum['subdive_features'],[1,0])
            _,_,posts_fine_1,_ = hhmm_FV_uncorr.fwd_bwd(datum['subdive_features'],[1,1])
            for i,(post_fine_0,post_fine_1) in enumerate(zip(posts_fine_0.T,posts_fine_1.T)):
                p = posts_crude.T[dive_num,0]*post_fine_0[1] + \
                    posts_crude.T[dive_num,1]*post_fine_1[1]
                data[dive_num]['subdive_features'][i]['hhmm_FV_uncorr_subdive'] = p
                data[dive_num]['subdive_features'][i]['hhmm_FV_uncorr_dive'] = posts_crude.T[dive_num,1]

        # get confusion matrix
        CM_coarse = np.zeros((2,2))
        CM_fine = [np.zeros((2,2)),np.zeros((2,2))]
        for dive in data:
            dive_type = dive['dive_type']
            CM_coarse[dive_type,1] += dive['subdive_features'][0]['hhmm_FV_uncorr_dive']
            CM_coarse[dive_type,0] += max(0,1.0-dive['subdive_features'][0]['hhmm_FV_uncorr_dive'])
            for seg in dive['subdive_features']:
                subdive_type = seg['subdive_type']
                CM_fine[dive_type][subdive_type,1] += seg['hhmm_FV_uncorr_subdive']
                CM_fine[dive_type][subdive_type,0] += max(0,1.0-seg['hhmm_FV_uncorr_subdive'])

        hhmm_FV_uncorr.CM = [CM_coarse,CM_fine]

        hhmm_FV_uncorr.save('../Params/hhmm_FV_uncorr_%d_%d'%(dataset_num,rand_seed))


    ### CarHHMM, no Z2 ###
    if model == 'CarHHMM1':
        print('')
        print('STARTING CarHHMM minux Z2')
        print('')
        pars = Parameters.Parameters()
        pars.K = [2,2]
        pars.features = [{'dive_duration':{'corr':False,'f':'gamma'}},
                         {'A':{'corr':True,'f':'normal'}}]

        hhmm_V = HHMM.HHMM(pars,data_V)
        hhmm_V.theta = hhmm_V_theta
        hhmm_V.eta = hhmm_V_eta
        hhmm_V.true_theta = deepcopy(hhmm_V_theta)
        hhmm_V.true_eta = deepcopy(hhmm_V_eta)

        hhmm_V.train_DM(data_V)
        hhmm_V.get_SEs(data_V,h)
        data,data_V,data_FV = create_data()

        # crude posterior
        _,_,posts_crude,_ = hhmm_V.fwd_bwd(data_V,[0])

        # fine posterior
        for dive_num,datum in enumerate(data_V):
            data[dive_num]['hhmm_v'] = posts_crude.T[dive_num,1]
            _,_,posts_fine_0,_ = hhmm_V.fwd_bwd(datum['subdive_features'],[1,0])
            _,_,posts_fine_1,_ = hhmm_V.fwd_bwd(datum['subdive_features'],[1,1])
            for i,(post_fine_0,post_fine_1) in enumerate(zip(posts_fine_0.T,posts_fine_1.T)):
                p = posts_crude.T[dive_num,0]*post_fine_0[1] + \
                    posts_crude.T[dive_num,1]*post_fine_1[1]
                if i%2 == 0:
                    p0 = p
                else:
                    data[dive_num]['subdive_features'][int((i-1)/2)]['hhmm_V_subdive'] = (p*p0)/(p*p0+(1.-p)*(1.-p0))
                    data[dive_num]['subdive_features'][int((i-1)/2)]['hhmm_V_dive'] = posts_crude.T[dive_num,1]

        # get confustion matrix
        CM_coarse = np.zeros((2,2))
        CM_fine = [np.zeros((2,2)),np.zeros((2,2))]
        for dive in data:
            dive_type = dive['dive_type']
            CM_coarse[dive_type,1] += dive['subdive_features'][0]['hhmm_V_dive']
            CM_coarse[dive_type,0] += max(0,1.0-dive['subdive_features'][0]['hhmm_V_dive'])
            for seg in dive['subdive_features']:
                subdive_type = seg['subdive_type']
                CM_fine[dive_type][subdive_type,1] += seg['hhmm_V_subdive']
                CM_fine[dive_type][subdive_type,0] += max(0,1.0-seg['hhmm_V_subdive'])

        hhmm_V.CM = [CM_coarse,CM_fine]
        hhmm_V.save('../Params/hhmm_V_%d_%d'%(dataset_num,rand_seed))


    ### CarHHMM ###
    if model == 'CarHHMM2':
        print('')
        print('STARTING CarHHMM')
        print('')
        pars = Parameters.Parameters()
        pars.K = [2,2]
        pars.features = [{'dive_duration':{'corr':False,'f':'gamma'}},
                         {'FoVeDBA':{'corr':False,'f':'gamma'},
                          'A':{'corr':True,'f':'normal'}}]

        hhmm_FV = HHMM.HHMM(pars,data_FV)
        hhmm_FV.theta = hhmm_FV_theta
        hhmm_FV.eta = hhmm_FV_eta
        hhmm_FV.true_theta = deepcopy(hhmm_FV_theta)
        hhmm_FV.true_eta = deepcopy(hhmm_FV_eta)

        hhmm_FV.train_DM(data_FV,max_iters=5)
        hhmm_FV.get_SEs(data_FV,h)
        data,data_V,data_FV = create_data()

        # crude posterior
        _,_,posts_crude,_ = hhmm_FV.fwd_bwd(data_FV,[0])

        # fine posterior
        for dive_num,datum in enumerate(data_FV):
            _,_,posts_fine_0,_ = hhmm_FV.fwd_bwd(datum['subdive_features'],[1,0])
            _,_,posts_fine_1,_ = hhmm_FV.fwd_bwd(datum['subdive_features'],[1,1])
            for i,(post_fine_0,post_fine_1) in enumerate(zip(posts_fine_0.T,posts_fine_1.T)):
                p = posts_crude.T[dive_num,0]*post_fine_0[1] + \
                    posts_crude.T[dive_num,1]*post_fine_1[1]
                data[dive_num]['subdive_features'][i]['hhmm_FV_subdive'] = p
                data[dive_num]['subdive_features'][i]['hhmm_FV_dive'] = posts_crude.T[dive_num,1]

        # get confusion matrix
        CM_coarse = np.zeros((2,2))
        CM_fine = [np.zeros((2,2)),np.zeros((2,2))]
        for dive in data:
            dive_type = dive['dive_type']
            CM_coarse[dive_type,1] += dive['subdive_features'][0]['hhmm_FV_dive']
            CM_coarse[dive_type,0] += max(0,1.0-dive['subdive_features'][0]['hhmm_FV_dive'])
            for seg in dive['subdive_features']:
                subdive_type = seg['subdive_type']
                CM_fine[dive_type][subdive_type,1] += seg['hhmm_FV_subdive']
                CM_fine[dive_type][subdive_type,0] += max(0,1.0-seg['hhmm_FV_subdive'])

        hhmm_FV.CM = [CM_coarse,CM_fine]

        # save hmms
        hhmm_FV.save('../Params/hhmm_FV_%d_%d'%(dataset_num,rand_seed))
