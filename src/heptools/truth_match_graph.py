# SwiftSimGNN
# ak -> networkx

"""
Important resource: https://stackoverflow.com/questions/70452465/how-to-load-in-graph-from-networkx-into-pytorch-geometric-and-set-node-features
"""

import awkward as ak
import networkx as nx
import uproot
import functools
import pyarrow as pa
import pandas as pd
import tqdm

def friend_from_ROOT(inputs,list_of_branches):

    list_of_branches.append("eventNumber")

    dfs = []

    for file in inputs:
        filename = file[0]
        treename = file[1]

        # Extract the eventNumber TBranch as a pandas DataFrame
        # Require this as separate for index?
        with uproot.open(f"{filename}:{treename}") as events:
            eventNumbers = events.arrays(["eventNumber"],library="pd")

        # Seems to work without line below so drop?
        eN = eventNumbers.reset_index(level=1, drop=True)

        # Extract the full desired Branches as awkward array, and parse
        with uproot.open(f"{filename}:{treename}") as events:
            outputs = events.arrays(list_of_branches,library="ak")
        arrow_array = ak.to_arrow(outputs)
        df = arrow_array.to_pandas()
        dff = df.to_frame()

        # Merge and set eventNumber as index
        merged = pd.concat([eN,dff],axis=1)
        final = merged.set_index(["eventNumber"])
        dfs.append(final)

    df_out = functools.reduce(lambda left, right: pd.merge(left, right, how="inner",on='eventNumber'), dfs)
    column_names = {df_out.columns[og]:new for og,new in enumerate([t[1] for t in inputs])}
    df_final = df_out.rename(columns=column_names)#, inplace=True)

    return df_final


def build_jet_event_graph(event):
    # event is a dictionary
    jet_events = []
    for i in range(len(event["jet_pt"])):
        ji = (1,event["jet_pt"][i],event["jet_eta"][i],event["jet_phi"][i],event["jet_e"][i])
        jet_events.append(ji)
    
    # Build graph
    G = nx.Graph()
    for i,jet in enumerate(jet_events):
        G.add_node(i,v=jet)

    GC = nx.complete_graph(G)

    return GC

def root_to_matched_graph(inputs,list_of_branches):

    #  Extract TTrees and turn into dataframe of linked events
    df_out = friend_from_ROOT(inputs,list_of_branches)

    tree_graph_dataframe = {}

    for index, row in tqdm(df_out.iterrows()):

        reco_event = row["reco"]
        pL_event   = row["particleLevel"]

        GC_reco = build_jet_event_graph(reco_event)
        GC_pL   = build_jet_event_graph(pL_event)

        tree_graph_dataframe[index] = {"particleLevel": GC_pL, "reco": GC_reco}

    Graph_DF = pd.DataFrame.from_dict(tree_graph_dataframe)

    return Graph_DF