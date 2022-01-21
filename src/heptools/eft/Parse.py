import pandas as pd
import json

"""
Combined parsers into single file using classes as namepsace - maybe bad Python practice?
TXT_Parse is untested here
"""


class JSON_Parse:

    @staticmethod
    def parse_limits(DL):
        return [(float(a),float(b)) for a,b in zip(DL[0],DL[1])]

    @classmethod
    def parse_from_json(self,filename):

        df = pd.io.json.read_json(filename).T

        df2 = df.drop(columns='Global mode')

        for k,v in df2.iterrows():
            v['1sigma Limits'] = self.parse_limits(v['1sigma Limits']) 
            v['2sigma Limits'] = self.parse_limits(v['2sigma Limits'])

        dft = pd.concat({"Bounds":df2.T}).T

        dft['Global mode'] = df['Global mode']

        return dft

    @classmethod
    def auto_parse(self,base_path,filename):
        dic_of_directories = {"Linear+Quadratic (Marg.)"     : base_path                + "/limits.json",
                                "Linear (Marg.)"             : base_path+"_linear"      + "/limits.json",
                                "Linear+Quadratic (Indp.)"   : base_path+"_independent" + "/limits.json"}

        d= {k:parse_from_json(fn) for k,fn in dic_of_directories.items()}
        
        return pd.concat(d, axis=1)





class TXT_Parse:

    @staticmethod
    def compute_mid_point(lower,upper): 
        return 0.5*(upper - lower)

    @staticmethod
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

    @staticmethod
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


    @classmethod
    def parse_general(self,bounds_filename,statistics_dir):

        """
        parse_general(...) extracts bounds from the given 1dlimits.txt file
        It also loops over all the files in the given directory with name containing 'statistics' and parse the global modes from these
        """

        overall_bounds_data = {}
        new_dict = {}
        if isinstance(bounds_filename,str):
            bounds_data       =  self.extract_bounds_from_txt(bounds_filename)
        if isinstance(bounds_filename,dict):
            # print("doing")
            for label,filename in bounds_filename.items():
                overall_bounds_data[label] = self.extract_bounds_from_txt(filename)

        df=pd.concat({"Bounds":pd.DataFrame.from_dict(overall_bounds_data).T}).T#,names=["FirstLevel"]).T

        import os 
        stats_files2parse = [file for file in os.listdir(statistics_dir) if "statistics" in file]
    
        # Loop over the statistics files
        global_mode_dict = {}
        for file in stats_files2parse:
            indiv_global_mode_data = self.extract_global_mode(statistics_dir +"/"+file)
            global_mode_dict = {**global_mode_dict,**indiv_global_mode_data}

        df["Global Mode"] = pd.Series(global_mode_dict)
        return df


    @classmethod
    def loop_dirs(self,dic_of_directories):
        # Eventually just define the base directory and everything is done from there
        d = {}
        for (k,v) in dic_of_directories.items():
            df = self.parse_general(bounds_filename={"68":v+"/Numerics/1dlimits.txt","95":v+"/Numerics/1dlimits_2sigma.txt"},statistics_dir=v+"/Numerics/")
            d[k] = df
        return pd.concat(d, axis=1)


    @classmethod
    def auto_parse(self,base_path):
        dic_of_directories = {"Linear+Quadratic (Marg.)"     : base_path,
                                "Linear (Marg.)"             : base_path+"_linear",
                                "Linear+Quadratic (Indp.)"   : base_path+"_independent"}
        
        output_df = self.loop_dirs(dic_of_directories)

        return output_df

