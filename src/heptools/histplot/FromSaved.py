
import os
from ROOT import TFile,TH1D

from heptools.histplot.PyHist_Class import Histogram_Wrapper as Hist
from heptools.histplot.Plotting import Ratio_Plot_ROOT



def ratio_plot_from_ROOT(list_of_ROOT_hists,divisor,**kwargs):

	dic_of_pyhists = {}
	for h in list_of_ROOT_hists:
		dic_of_pyhists[h.GetName()] = Hist(h,name=h.GetName(),colour="black" ,legend_entry=h.GetName())
	
	py_divisor = dic_of_pyhists[divisor.GetName()]
	p = Ratio_Plot_ROOT("A Plot",list_of_histograms=[x_nom,x_up,x_down],divisor=x_nom,normalise=False)
	p.Initialise_Plot_Design("ATLAS")
	plt,ax,rax = p.Make_Ratio_Plot("basic-line")

	
ratio_plot_from_ROOT([1,2,3],1)



def PyHist_fromROOT(filename,histname,**kwargs):

	'''
	Function extracts ROOT TTree from associated ROOT file
	Could implement a sub-function here such that it can be used to loop over many files
	'''

	assert os.path.isfile(filename), "Cannot find this file: "+ filename
	file = TFile(filename,"READ")

	ROOT_hist = file.Get(histname)

	normalise = kwargs["normalise"] if "normalise" in kwargs else True

	if normalise:
		ROOT_hist.Scale(1/ROOT_hist.Integral())

	ROOT_hist.SetDirectory(0)

	return ROOT_hist









