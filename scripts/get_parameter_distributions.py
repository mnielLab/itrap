#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import scipy.stats as st
from scipy.stats._continuous_distns import _distn_names
import warnings

# https://docs.scipy.org/doc/scipy/reference/tutorial/stats.html
# https://stackoverflow.com/questions/6620471/fitting-empirical-distribution-to-theoretical-ones-with-scipy-python?lq=1
# https://stackoverflow.com/questions/37487830/how-to-find-probability-distribution-and-parameters-for-real-data-python-3
#########################################################################################################
#                                               Functions                                               #
#########################################################################################################
class norm2_gen(st.rv_continuous):
    def _argcheck(self, *args):
        return True

    def _pdf(self, x, m, s, w, delta):
        phi = 0.5 + np.arctan(w)/np.pi
        return np.exp(-(x-m+delta/2)**2 / (2. * s**2)) / np.sqrt(2. * np.pi * s**2) * phi + \
               np.exp(-(x-m-delta/2)**2 / (2. * s**2)) / np.sqrt(2. * np.pi * s**2) * (1 - phi)
    
# Create models from data
def best_fit_distribution(data, bins=200, ax=None):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Best holders
    best_distributions = []
    
    distribution_collection = [d for d in _distn_names if not d in ['levy_stable', 'studentized_range']] #+ ['norm2']
    #distribution_collection = [d for d in dir(st) if isinstance(getattr(st, d), st.rv_continuous)]

    # Estimate distribution parameters from data
    for ii, distribution_name in enumerate(distribution_collection):

        print("{:>3} / {:<3}: {}".format( ii+1, len(distribution_collection), distribution_name ))

        #distribution = getattr(st, distribution_name)
        try:
            distribution = getattr(st, distribution_name)
        except AttributeError:
            distribution = norm2_gen(name=distribution_name)

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                
                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]
                
                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))
                
                # if axis pass in add to plot
                try:
                    if ax:
                        pd.Series(pdf, x).plot(ax=ax)
                    end
                except Exception:
                    pass

                # identify if this distribution is better
                #print(distribution_name, sse)
                best_distributions.append((distribution, params, sse))
        
        except Exception:
            pass

    return sorted(best_distributions, key=lambda x:x[2])

def make_pdf(dist, params, data, size=10000):
    """Generate distributions's Probability Distribution Function """
    
    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    #start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    #end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)
    start = 0
    end = np.ceil(data.max())+1

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf

def plot_data(variable):
    data = df[variable].dropna()
    
    if data.fillna(0).apply(float.is_integer).all():
        bins = int(data.max()+1)
    else:
        bins = int(data.max()+1) * 10
    
    # Plot for comparison
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(20,6))
    data.plot(ax=ax1, kind='hist', bins=bins, density=True, alpha=0.5, #bins=50
              color=list(matplotlib.rcParams['axes.prop_cycle'])[1]['color'])

    # Save plot limits
    dataYLim = ax1.get_ylim()
    dataXLim = ax1.get_xlim()

    # Find best fit distribution
    best_distibutions = best_fit_distribution(data, bins, ax1)
    best_dist = best_distibutions[0]

    # Update plots
    ax1.set_ylim(dataYLim)
    ax1.set_title(u'All Fitted Distributions')
    ax1.set_xlabel(u'%s' %variable.replace('_',' '))
    ax1.set_ylabel('Frequency')

    # Make PDF with best params 
    pdf = make_pdf(best_dist[0], best_dist[1], data)

    # Display
    pdf.plot(ax=ax2, lw=2, label='PDF', legend=True)
    data.plot(kind='hist', bins=50, density=True, alpha=0.5, label='Data', legend=True, ax=ax2)

    param_names = ((best_dist[0].shapes + ', loc, scale').split(', ')
                   if best_dist[0].shapes else ['loc', 'scale'])
    param_str = ', '.join(['{}={:0.2f}'.format(k,v) for k,v in zip(param_names, best_dist[1])])
    dist_str = '{}({})'.format(best_dist[0].name, param_str)

    ax2.set_title(u'Best fit distribution \n%s' %dist_str)
    ax2.set_xlabel(u'%s' %variable.replace('_',' '))
    ax2.set_ylabel('Frequency')
    ax2.set_xlim(dataXLim) #dataXLim
    ax2.set_ylim(dataYLim)
    plt.suptitle('%s' %variable.replace('_',' '))
    for f in PLOT:
        plt.savefig(f %variable, bbox_inches='tight')
    
    return best_dist[0].name, best_dist[1]

parameters = ['umi_count_mhc','delta_umi_mhc','umi_count_mhc_rel',
              'umi_count_cd8','delta_umi_cd8',
              'umi_count_TRA','delta_umi_TRA',
              'umi_count_TRB','delta_umi_TRB']

#########################################################################################################
#                                                 Input                                                 #
#########################################################################################################

INPUT = snakemake.input.data #
OUTPUT = snakemake.output.prms
PLOT = snakemake.params.plot

#########################################################################################################
#                                                 Load                                                  #
#########################################################################################################
df = pd.read_csv(INPUT)
#df.fillna(dict.fromkeys(parameters, 0), inplace=True)

#########################################################################################################
#                                                 Main                                                  #
#########################################################################################################
dct = dict()
for parameter in parameters:
    print(parameter)
    name, args = plot_data(parameter)
    dct[parameter] = [name, args]

pd.DataFrame.from_dict(dct, orient='index', columns=['dist_name','dist_args']).to_csv(OUTPUT)

