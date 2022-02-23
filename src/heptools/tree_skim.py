# Python3.6 ROOT macro for splitting TTree

# from ctypes.wintypes import HHOOK
# from ROOT import TFile, TTree

from ROOT import TTree, TFile, vector, gROOT, TLorentzVector


import sys

# sys.argv[]

# _file0 = TFile("mc15_13TeV.410011.PwPyEG_P2012_singletop_tchan_lept_top.1lep_raw.root")

# tree = _file0.Get("mini")

def ResetBranches(dictionary):

	for key in dictionary:
		dictionary[key].resize(0)


def BookBranch(name, suffix, tree):

    """ Function to book branches to an output TTree """

    branch_type = 'null'
    if "F" in suffix:
        branch_type = 'float'
    elif "I" in suffix:
        branch_type = 'int'
    else:
        print("WARNING: Type not recognised",suffix)

    branch = vector(branch_type)(20)
    tree.Branch(name, branch)
    return branch

def skim_tree(input_tree,branches2keep):

    # Define new tree name
    new_tree_name = input_tree.GetName() + "_skimmed"
    # Open new output file
    output_file = TFile(new_tree_name+".root", "RECREATE")

    # Output TTree
    output_tree   = TTree(new_tree_name,   new_tree_name)

    # Turn off all irrelevant branches, and turn back on the pertinent ones
    input_tree.SetBranchStatus("*",0)
    for i in branches2keep:
        input_tree.SetBranchStatus(i, 1)
        
    # Store branches here
    branch_dict = {} 
    for branch in branches2keep:
        branch_dict[branch[0]] = BookBranch(branch[0], branch[1], output_tree)        

    total_events = input_tree.GetEntries()

    for counter, event in enumerate(input_tree):
        
        if counter % 1000 == 0: 
            print("Event",counter,"/",total_events)
            
        ResetBranches(branch_dict) #Don't tell me what to do!
        
        input_tree.GetEntry(counter) #this
        
        for branch in branches2keep:
            branch_dict["weight_mc"].push_back(getattr(event,branch[0]))

