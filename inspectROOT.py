from ROOT import TFile

import sys

path = sys.argv[1]

tree_names = []
 
tfile = TFile.Open(path)

for key in tfile.GetListOfKeys():
	print(key)