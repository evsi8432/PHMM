import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

import sys

from scipy.stats import gamma
from scipy.stats import norm
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde
from scipy.special import expit
from scipy.special import logsumexp
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from scipy.signal import convolve
from scipy.interpolate import interp1d

import divebomb

import time
import pickle

import Preprocessor
import Parameters
import HHMM
import Visualisor

np.random.seed(0)

model = str(sys.argv[1])

# set parameters
pars = Parameters.Parameters()
pars.features = [{'dive_duration':{'corr':False,'f':'gamma'}},
                 {'Ahat_low':{'thresh':5,'corr':False,'f':'gamma'},
                  'Ax':{'corr':True,'f':'normal'},
                  'Ay':{'corr':True,'f':'normal'},
                  'Az':{'corr':True,'f':'normal'}}]
pars.K = [2,3]
pars.share_fine_states = True

if model == 'CarHHM':
    pars.K = [1,3]
elif model == 'HHMM':
    pars.features[1]['Ax']['corr'] = False
    pars.features[1]['Ay']['corr'] = False
    pars.features[1]['Az']['corr'] = False
elif model == 'CarHHMM1':
    pars.features[1] = {'Ax':{'corr':True,'f':'normal'},
                        'Ay':{'corr':True,'f':'normal'},
                        'Az':{'corr':True,'f':'normal'}}

# do preprocessing
prep = Preprocessor.Preprocessor(pars)
data_outfile = '../Data/Data_%s_k_%s_%s_dives_same_fine_states' % (model,pars.K[0],pars.K[1])

print('loading data')
df = prep.load_data(pars.cvc_file,pars.csv_file,pars.cvc_cols)
print('pruning cols')
df = prep.prune_cols(df)
print('pruning times')
df = prep.prune_times(df,pars.stime,pars.etime,pars.drop_times)
print('fixing pressure')
df = prep.fix_pressure(df)
print('finding Vz')
df = prep.find_Vz(df)
print('smoothing cols')
df = prep.smooth_columns(df,pars.smoother,pars.smooth_cols)
print('dividing into dives')
df,dive_df = prep.find_dives(df)
print('getting features')
data = prep.get_all_features(df,dive_df)

with open(data_outfile, 'wb') as f:
    pickle.dump(data, f)
