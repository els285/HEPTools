import matplotlib.pyplot as plt
import numpy as np
from ROOT import TFile,TAxis,TH1,gROOT
import os
import numpy as np
import pickle

# TH1.SetDefaultSumw2()



class PyHist:
    """ Basic wrapper for ROOT histogram 
    Should contain no ROOT functionality, just a container for the information
    """

    def __init__(self,Name,Bin_Values,Bin_Errors,Bin_Centres,Bin_Edges,**kwargs):

        # Data
        self.Name        = Name
        self.Bin_Values  = Bin_Values  
        self.Bin_Errors  = Bin_Errors  
        self.Bin_Centres = Bin_Centres 
        self.Bin_Edges   = Bin_Edges   

         # Plotting information
        self.legend_entry       = kwargs["legend_entry"]     if "legend_entry"      in kwargs else self.Name
        self.design_dict = kwargs["design_dict"]



class Histogram_Wrapper:
    """
    Larger wrapper which contains PyHists, including methods for ROOT->numpy conversion
    """

    @staticmethod
    def get_bin_values(hist):
        return np.asarray([hist.GetBinContent(binn+1) for binn in range(0,hist.GetNbinsX())])

    @staticmethod
    def get_bin_errors(hist):
        return np.asarray([hist.GetBinError(binn+1) for binn in range(0,hist.GetNbinsX())])

    @staticmethod
    def get_bin_edges(hist):
        return np.asarray([hist.GetXaxis().GetBinUpEdge(binn) for binn in range(0,hist.GetXaxis().GetNbins()+1)])

    @staticmethod
    def get_bin_centres(hist):
        return np.asarray([hist.GetXaxis().GetBinCenter(binn+1) for binn in range(0,hist.GetNbinsX())])

    @staticmethod
    def Compute_Normalised(ROOT_hist):

        """ Takes a ROOT histogram and normalised it"""

        h1dN = ROOT_hist.Clone(ROOT_hist.GetName()+"_norm")
        h1dN.Scale(1/ROOT_hist.Integral())

        return h1dN

    @classmethod
    def Create_Wrapper(self,hist,name,**kwargs):

        Bin_Values  = self.get_bin_values(hist)
        Bin_Errors  = self.get_bin_errors(hist)
        Bin_Centres = self.get_bin_centres(hist)
        Bin_Edges   = self.get_bin_edges(hist)
        return PyHist(name,Bin_Values,Bin_Errors,Bin_Centres,Bin_Edges,**kwargs)

    def __init__(self,ROOT_hist,**kwargs):

        '''
        The Histogram_Wrapper object is for converting ROOT histograms into numpy objects,
            but the data members of the class should contain no ROOT objects.
        '''
        # Meta-data
        self.name               = kwargs["name"]             if "name"              in kwargs else ROOT_hist.GetName()
        self.observable_type    = kwargs["obs"]              if "obs"               in kwargs else None
        self.file_name          = kwargs["filename"]         if "filename"          in kwargs else None

        # Plotting information
        self.legend_entry       = kwargs["legend_entry"]     if "legend_entry"      in kwargs else self.name
        
        self.design_dict = {}
        self.design_dict["colour"]             = kwargs["colour"]           if "colour"            in kwargs else None
        self.design_dict["linewidth"]          = kwargs["linewidth"]        if "linewidth"         in kwargs else 2
        self.design_dict["line_style"]         = kwargs["line_style"]       if "line_style"        in kwargs else "-"
        self.design_dict["marker_style"]       = kwargs["marker_style"]     if "marker_style"      in kwargs else ""
        self.design_dict["opacity"]            = kwargs["opacity"]          if "opacity"           in kwargs else 0.4
        self.design_dict["error_opacity"]      = kwargs["error_opacity"]    if "error_opacity"     in kwargs else 0.2

        # ROOT histograms
        self.UnNorm_ROOT_hist = ROOT_hist


        # Unnormalised Hist wrapper
        self.UnNorm_PyWrap_Hist = self.Create_Wrapper(ROOT_hist,self.name+"_Unnormalised",legend_entry=self.legend_entry,design_dict=self.design_dict)

        # Normalising histograms should depend on if histogram is empty or not
        if ROOT_hist.Integral()!=0.0:
            norm_hist = self.Compute_Normalised(ROOT_hist)
            self.Norm_ROOT_hist   = norm_hist 
            # Normalised Hist wrapper
            self.Norm_PyWrap_Hist   = self.Create_Wrapper(norm_hist,self.name+"_Normalised",legend_entry=self.legend_entry,design_dict=self.design_dict)





