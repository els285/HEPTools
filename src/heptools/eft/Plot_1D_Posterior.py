# For plotting 1D posterior distribution (along with global model and another ancillary data)
from collections import namedtuple
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import boost_histogram as bh
import functools
import operator


def normalise_boost_histogram(boost_hist):

    area = functools.reduce(operator.mul, boost_hist.axes.widths)
    factor = np.sum(boost_hist.view())
    view = boost_hist.view() / (factor * area)

    for i,x in enumerate(view):  
        boost_hist[i]=x

    return boost_hist



class _1D_Posterior_Plot:

    bin_width = 0.1

    def make_binning(self,low,high):
        return np.concatenate([np.arange(low,high,self.bin_width)  ,  np.asarray([high])])

    def __init__(self,posterior,weights,limits,**kwargs):

        # Inputs are posterior data and limits
        self.posterior = posterior  # This should be a numpy array 
        self.weights   = weights
        self.limits    = limits.Bounds
        self.global_mode = limits['Global mode']
    
        # Use the maximum and minimum values in the posterior as bounds of the histogram
        mini=round(min(self.posterior),1)
        maxi=round(max(self.posterior),1)

        self.extrema = namedtuple('Extrema', ['mini', 'maxi'])(mini, maxi)

        # Use boost_histogram to generate the histogram object using the posterior data
        number_of_bins = 120
        hist = bh.Histogram(bh.axis.Regular(number_of_bins,mini,maxi))
        hist.fill(self.posterior,weight=self.weights)

        # Normalise the histogram
        hist_density = normalise_boost_histogram(hist)
        self.Primary_Hist =hist_density


        # The individual regions are defined as sliced histograms, from lower to upper bounds
        if len(self.limits["1sigma Limits"])!=0:
            self._1Sigma_Histograms_bh = [self.Primary_Hist[ bh.loc(x1.lower) : bh.loc(x1.upper)] for x1 in self.limits["1sigma Limits"]]

        if len(self.limits["2sigma Limits"])!=0:
            self._2Sigma_Histograms_bh = [self.Primary_Hist[ bh.loc(x2.lower) : bh.loc(x2.upper)] for x2 in self.limits["2sigma Limits"]]

        if "3sigma Limits" in self.limits.index and len(self.limits["3sigma Limits"])!=0:
            self._3Sigma_Histograms_bh = [self.Primary_Hist[ bh.loc(x3.lower) : bh.loc(x3.upper)] for x3 in self.limits["3sigma Limits"]]
        else:
            self._3Sigma_Histograms_bh = [ self.Primary_Hist[bh.loc(mini) : bh.loc(maxi)]]
        
        self.Region_Histograms = [self._3Sigma_Histograms_bh , self._2Sigma_Histograms_bh , self._1Sigma_Histograms_bh]
          


    def plot_posterior(self):

        for color,hist_list in zip(["red","yellow","lawngreen"] , self.Region_Histograms):
            for hist in hist_list:
                hep.histplot(hist.to_numpy(),fill=True,color=color)

        # # Black outline
        hep.histplot(self.Primary_Hist.to_numpy(),fill=False,color="black",linewidth=0.4)

        self.fig = plt.gcf() 
        self.ax  = plt.gca()


        self.make_legend()


    def make_legend(self):
        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D


        _1sigma_patch = mpatches.Patch(color='lawngreen', label=r'$1 \sigma $')
        _2sigma_patch = mpatches.Patch(color='yellow', label=r'$2 \sigma $')
        _3sigma_patch = mpatches.Patch(color='red', label=r'$3 \sigma $')

        global_line = Line2D([], [], color='black',  linestyle='dashed', label='Global mode')
        
        handles = [_1sigma_patch,_2sigma_patch,_3sigma_patch]

        if self.include_global_mode: handles.append(global_line)
        plt.legend(handles=handles)


    def add_global_mode(self):
        posterior_apex = max(self.Primary_Hist.to_numpy()[0])
        self.ax.vlines(self.global_mode,0,posterior_apex*1.08,color="black",linestyle='dashed',lw=0.85)
        
        from matplotlib.lines import Line2D

    
    def make_plot(self,include_global_mode,include_prior):

        self.include_global_mode = include_global_mode
        self.include_prior       = include_prior

        self.plot_posterior()
        if self.include_global_mode: self.add_global_mode()
        # if self.include_prior:       self.add_prior()
        self.make_legend()




class _1D_LHC_Posterior_Plot(_1D_Posterior_Plot):


    allowed_experiment_styles = ["ATLAS","ALICE","LHCb2","CMS","ATLASTex"]

    def initialise_LHC_plot(self,experiment):

        assert any([experiment == x for x in self.allowed_experiment_styles]), "Experiment style not defined"

        if   experiment=="ATLAS": plt.style.use(hep.style.ATLAS) # This is the correct syntax for Pytohn3.6.9 version (mplhep 0.2.8)
        elif experiment=="CMS"  : plt.style.use(hep.style.CMS)
        elif experiment=="ALICE": plt.style.use(hep.style.ALICE)
        elif experiment=="LHCb2": plt.style.use(hep.style.LHCb2)
        elif experiment=="ATLASTex": 
            print("This is not supported because do not have a working TexLive distribution")
            # plt.style.use(hep.style.ATLASTex) # T

    def include_metadata(self,text,data,lumi,year):
        hep.atlas.label(text, data=data, lumi=lumi, year=year)

    def additional_label(self,text,**kwargs):
        plt.text(0.05,0.785,text,transform=self.ax.transAxes)


    def __init__(self,posterior,weights,limits,experiment,**kwargs,):
        
        super().__init__(posterior,weights,limits,**kwargs)
        self.experiment = experiment

        self.initialise_LHC_plot(experiment)


        

        


