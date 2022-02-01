from heptools.skimming import generic_tree_skim
from ROOT import TFile


def main():

    _file0 = TFile("~/mc15_13TeV.410011.PwPyEG_P2012_singletop_tchan_lept_top.1lep_raw.root")

    tree = _file0.Get("mini")

    # branches2keep = [["eventNumber","I"],["lep_type","VF"],["lep_pt","VF"]]

    branches2keep = ["eventNumber","lep_type","lep_pt"]

    generic_tree_skim(tree,branches2keep,cut_off=200000,file_name="hi.root",tree_name="parton_tree")
    
main()