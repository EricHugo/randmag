#!/usr/bin/env python

"""
Fits distributions onto given dataset
"""

import operator
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy import stats as st

class FitDist():
    def __init__(self, dataset):
        # put list of data into pandas series
        self.dataset = pd.Series(tuple(dataset))
        print(self.dataset.value_counts())
    
    def histogram(self, distfit):
        """Produces histogram with desired fit and number of bins"""
        sns.distplot(self.dataset, fit=distfit)
        plt.savefig('distplot.png', transparent=True)
        #plt.show()
        return

    def find_best_dist(self):
        """Determines best fit for a number of distributions based
        on maxiumum likelihood"""
        distributions = [st.norm, st.gamma, st.beta]
        mles = {}
        for distribution in distributions:
            pars = distribution.fit(self.dataset)
            mle = distribution.nnlf(pars, self.dataset)
            mles[distribution] = mle
        #print(mles)
        best_fit = max(mles.items(), key=operator.itemgetter(1))[0]
        return best_fit

    def placholder(self, nan):
        pass
