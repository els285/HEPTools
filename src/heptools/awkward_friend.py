import awkward as ak
import awkward_pandas as akpd
import numpy as np
import pandas as pd
import uproot
from functools import reduce


def friend_from_ROOT(inputs,list_of_branches):

    """
    Build a friended pandas dataframe of awkward arrays
    Uses new uproot5 + awkward-pandas
    """
 
    list_of_branches.append("eventNumber")

    input_series = []

    for file in inputs:
        filename = file[0]
        treename = file[1]
        
        with uproot.open(f"{filename}:{treename}") as events:
            outputs = events.arrays(list_of_branches,library="ak")

        # Transform to pandas Series
        s = akpd.from_awkward(outputs)
        s.name = treename

        # Extract eventNumbers
        eventNumbers_array = outputs["eventNumber"].to_numpy()
        eventNumbers = np.array(list(map(lambda x: x[0], eventNumbers_array)))

        s_out = s.set_axis(eventNumbers)

        input_series.append(s_out)

    common_keys = reduce(np.intersect1d,[so.index for so in input_series])

    dfs = []
    for s in input_series:
        s2 = s[s.index.isin(common_keys)]
        s3 = s2.to_frame()
        s3.index.name = "eventNumber"
        dfs.append(s3)

    df_out = reduce(lambda left, right: pd.merge(left, right, how="inner",on='eventNumber'), dfs)

    return df_out