from tqdm import tqdm
import uproot
import numpy as np

# Edit this data
    
root_file       = "/Users/ethansimpson/delphes_output_0706_ttbar_madlad.root"
TTree_name      = "Delphes"
output_filename = "delphes_matched1.root"
output_treename = "tree" 
    
# Useful functions
nested_hstack = lambda L:  [q for p in L for q in p]
DeltaR        = lambda X,Y: (X**2 + Y**2)**(0.5)

def matching(list_of_track_eta,list_of_track_phi,jet_eta,jet_phi):

    """
    Returns the matched indices of the tracks corresponding to the jet in question
    """
    
    eta_diff = list_of_track_eta - jet_eta*np.ones(len(list_of_track_eta))
    phi_diff = list_of_track_phi - jet_phi*np.ones(len(list_of_track_phi))
    DR_array = DeltaR(eta_diff,phi_diff)
    return np.argwhere(DR_array<0.4)[:,0]


def process(jet_tree,track_tree):
    output_dic = {}
    Nevents = len(jet_tree)
    for eN in tqdm(range(Nevents)):
        inner_dic = {}
        Njets = len(jet_tree[eN]["Jet.Eta"])
        track_Eta_arr = track_tree["Track.Eta"][eN]
        track_Phi_arr = track_tree["Track.Phi"][eN]
        for jet_index in range(Njets):
            jet_eta = jet_tree[eN]["Jet.Eta"][jet_index]
            jet_phi = jet_tree[eN]["Jet.Phi"][jet_index]
            matched_indices = matching(list_of_track_eta=track_Eta_arr,
                    list_of_track_phi=track_Phi_arr,
                    jet_eta=jet_eta,jet_phi=jet_phi)
            inner_dic[jet_index] = {"jet_pt"    :   jet_tree[eN]["Jet.PT"][jet_index],
                                    "jet_eta"   :   jet_tree[eN]["Jet.Eta"][jet_index],
                                    "jet_phi"   :   jet_tree[eN]["Jet.Phi"][jet_index],
                                    "jet_mass"  :   jet_tree[eN]["Jet.Mass"][jet_index],
                                    "track_pt"  :   track_tree[eN]["Track.PT"][matched_indices],
                                    "track_eta" :   track_tree[eN]["Track.Eta"][matched_indices],
                                    "track_phi" :   track_tree[eN]["Track.Phi"][matched_indices]
                                    }
        output_dic[eN] = inner_dic   
    return output_dic

    
  
def restructure_data(output_dic):

    rearranged = {  "Njets"     : [],
                    "jet_pt"    : [],
                    "jet_eta"   : [],
                    "jet_phi"   : [],
                    "jet_mass"  : [],
                    "track_pt"  : [],
                    "track_eta" : [],
                    "track_phi" : [] }
    for event in output_dic.values():
        jet_pt      = []
        jet_eta     = []
        jet_phi     = []
        jet_mass    = []
        track_pt    = []
        track_eta   = []
        track_phi   = []
        rearranged["Njets"].append(len(event))
        for jet in event.values():
            jet_pt.append(jet["jet_pt"])
            jet_eta.append(jet["jet_eta"])
            jet_phi.append(jet["jet_phi"])
            jet_mass.append(jet["jet_mass"])
            track_pt.append(jet["track_pt"])
            track_eta.append(jet["track_eta"])
            track_phi.append(jet["track_phi"])
        rearranged["jet_pt"].append(jet_pt)
        rearranged["jet_eta"].append(jet_eta)
        rearranged["jet_phi"].append(jet_phi)
        rearranged["jet_mass"].append(jet_mass)
        rearranged["track_pt"].append(track_pt)
        rearranged["track_eta"].append(track_eta)
        rearranged["track_phi"].append(track_phi)
        
        return rearranged
    
def write_output(rearranged,output_filename,output_treename):
    eventNumbers_array = [i*np.ones(M) for i,M in enumerate(rearranged["Njets"])]
    with uproot.recreate(f"{output_filename}") as file:
        file[f"{output_treename}"] = {"eventNumber": np.hstack(eventNumbers_array),
                        "jet_pt"   : np.hstack(rearranged["jet_pt"]),
                        "jet_eta"   : np.hstack(rearranged["jet_eta"]),
                        "jet_phi"   : np.hstack(rearranged["jet_phi"]),
                        "jet_mass"  : np.hstack(rearranged["jet_mass"]),
                        "track_pt"  : nested_hstack(rearranged["track_pt"]),
                        "track_eta" : nested_hstack(rearranged["track_eta"]),
                        "track_phi" : nested_hstack(rearranged["track_phi"])}  
    print("Output written")      


def main():

    
    loaded_tree = uproot.open(f"{root_file}:{TTree_name}")
    jet_tree   = loaded_tree.arrays(["Jet.PT","Jet.Eta","Jet.Phi","Jet.Mass"],library="ak")
    track_tree = loaded_tree.arrays(["Track.PT","Track.Eta","Track.Phi"],library="ak")
    
    matched_data_dic = process(jet_tree,track_tree)
    output           = restructure_data(matched_data_dic)
    write_output(output,output_filename,output_treename)
      
      
main()