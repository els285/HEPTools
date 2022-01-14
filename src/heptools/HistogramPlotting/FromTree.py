import os
from ROOT import TFile,TH1D



def from_saved(filename,histname,**kwargs):

	'''
	Function extracts ROOT TTree from associated ROOT file
	Could implement a sub-function here such that it can be used to loop over many files
	'''

	assert os.path.isfile(filename), "Cannot find this file: "+ filename
	file = TFile(filename,"READ")

	ROOT_hist = file.Get(histname)

	normalise = kwargs["normalise"] if "normalise" in kwargs else False

	if normalise:
		ROOT_hist.Scale(1/ROOT_hist.Integral())

	ROOT_hist.SetDirectory(0)

	return ROOT_hist




def project_hist_unNormalised(filename,treename,binning,branchname):

	"""
	Simple command to pull unnormalised histograms from ROOT TTrees
	"""

	assert os.path.isfile(filename), "Cannot find this file: "+ filename
	file = TFile(filename,"READ")
	assert hasattr(file,treename), "In file "+filename+", TTree not found:" + treename
	tree = file.Get(treename)

	ROOT_hist = TH1D(branchname+"_1D","",len(binning)-1,binning)
	tree.Project(ROOT_hist.GetName(),branchname)

	# ROOT_hist.Scale(1/ROOT_hist.Integral())
	ROOT_hist.SetDirectory(0)

	return ROOT_hist









