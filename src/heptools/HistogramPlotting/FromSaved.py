
import os
from ROOT import TFile,TH1D
# from PyHep_Hist_no_rootnumpy import Histogram



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









