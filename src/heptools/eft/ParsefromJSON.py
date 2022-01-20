import json
import pandas as pd
from collections import namedtuple


def parse_limits(DL):
    Pair_Bounds = namedtuple('Bounds_Tuple','lower upper')
    return [Pair_Bounds(float(a),float(b)) for a,b in zip(DL[0],DL[1])]


def parse_from_json(filename):

    df = pd.io.json.read_json(filename).T

    df2 = df.drop(columns='Global mode')

    for k,v in df2.iterrows():
        v['1sigma Limits'] = parse_limits(v['1sigma Limits']) 
        v['2sigma Limits'] = parse_limits(v['2sigma Limits'])

    dft = pd.concat({"Bounds":df2.T}).T

    dft['Global mode'] = df['Global mode']

    return dft


def auto_parse(base_path,filename):
    dic_of_directories = {"Linear+Quadratic (Marg.)"     : base_path                + "/limits.json",
                            "Linear (Marg.)"             : base_path+"_linear"      + "/limits.json",
                            "Linear+Quadratic (Indp.)"   : base_path+"_independent" + "/limits.json"}
    

    d= {k:parse_from_json(fn) for k,fn in dic_of_directories.items()}

    return pd.concat(d, axis=1)

# filename = "/home/ethan/EFTfitterSpinCorr.jl/nov15_eft_limits.json"

# print(parse_from_json(filename))


