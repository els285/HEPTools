from ROOT import TFile

_file0 = TFile("mini_skimmed.root")

tree = _file0.Get("mini_skimmed")


for event in tree:
    print(event.lep_type)
    input()