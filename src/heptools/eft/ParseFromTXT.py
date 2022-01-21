# For parsing a series of text files 

"""
Basic plotting format
cHQ1	[-5.2],[7.0]
cHQ3	[-8.4],[4.0]
cHt	[-4.6],[8.4]
"""

import re 
import pandas as pd
import matplotlib.pyplot as plt


def compute_mid_point(lower,upper): 
    return 0.5*(upper - lower)


def extract_bounds_from_txt(filename):

    data = {}
    file1 = open(filename, 'r')
    Lines = file1.readlines()

    for x in Lines:
        indiv_line = x.replace("\t",",").replace("\n","")
        t = indiv_line.split("[")
        wc=t[0].replace(",","")
        lower_bounds_str_list = t[1][:-1].replace("]","").split(",")
        upper_bounds_str_list = t[2].replace("]","").split(",")
        bounds = [(float(l),float(u)) for (l,u) in zip(lower_bounds_str_list,upper_bounds_str_list)]
        data[wc] = bounds

    return data


def extract_global_mode(filename):
    """
    The lines in the txt file are NamedTuples, which are a Python object. Surely we can push to some other filetype from that?
    """
    
    data = {}
    file1 = open(filename, 'r')
    Lines = file1.readlines()

    for x in Lines:
        l1 =x.split("[")[1].replace("\n","").split(")")
        for m in l1:
            l2 = m.strip("").replace("(","").split(",")
            for x in l2:
                if "parameter" in x:
                    wilson_coef = x.replace(" ","").split("=")[1].replace(":","")
                if "global_mode" in x:
                    global_mode = x.replace(" ","").split("=")[1]
            assert "wilson_coef" in locals(), "parameter term not found"
            assert "global_mode" in locals(),"global mode value not found"
            data[wilson_coef] = float(global_mode)
    
    return data

# def Parse(bounds_filename,statistics_filename):

#     bounds_data       =  extract_bounds_from_txt(bounds_filename)
#     global_mode_data  =  extract_global_mode(statistics_filename)

#     data = {}
#     for wc in bounds_data:
#         assert wc in global_mode_data.keys(),"Cannot find "+wc+" in global_mode_data"
#         data[wc] = {"Bounds":bounds_data[wc],"Global Mode":global_mode_data[wc]}

#     df = pd.DataFrame.from_dict(data).T
#     return df 


def parse_general(bounds_filename,statistics_dir):

    """
    parse_general(...) extracts bounds from the given 1dlimits.txt file
    It also loops over all the files in the given directory with name containing 'statistics' and parse the global modes from these
    """

    overall_bounds_data = {}
    new_dict = {}
    if isinstance(bounds_filename,str):
        bounds_data       =  extract_bounds_from_txt(bounds_filename)
    if isinstance(bounds_filename,dict):
        # print("doing")
        for label,filename in bounds_filename.items():
            overall_bounds_data[label] = extract_bounds_from_txt(filename)

    df=pd.concat({"Bounds":pd.DataFrame.from_dict(overall_bounds_data).T}).T#,names=["FirstLevel"]).T

    import os 
    stats_files2parse = [file for file in os.listdir(statistics_dir) if "statistics" in file]
   
    # Loop over the statistics files
    global_mode_dict = {}
    for file in stats_files2parse:
        indiv_global_mode_data = extract_global_mode(statistics_dir +"/"+file)
        global_mode_dict = {**global_mode_dict,**indiv_global_mode_data}


    df["Global Mode"] = pd.Series(global_mode_dict)
    return df
    # print(df)
    # input()

    # # df1 = df1.assign(e=pd.Series(np.random.randn(sLength)).values)

    # data = {}
    # for wc in bounds_data:
    #     assert wc in global_mode_dict.keys(),"Cannot find "+wc+" in global_mode_data"
    #     data[wc] = {"Bounds":bounds_data[wc],"Global Mode":global_mode_dict[wc]}

    # df = pd.DataFrame.from_dict(data).T
    # return df 


def loop_dirs(dic_of_directories):
    # Eventually just define the base directory and everything is done from there
    d = {}
    for (k,v) in dic_of_directories.items():
        df = parse_general(bounds_filename={"68":v+"/Numerics/1dlimits.txt","95":v+"/Numerics/1dlimits_2sigma.txt"},statistics_dir=v+"/Numerics/")
        d[k] = df
    return pd.concat(d, axis=1)


def auto_parse(base_path):
    dic_of_directories = {"Linear+Quadratic (Marg.)"     : base_path,
                            "Linear (Marg.)"            : base_path+"_linear",
                            "Linear+Quadratic (Indp.)"    : base_path+"_independent"}
    
    output_df = loop_dirs(dic_of_directories)

    return output_df

