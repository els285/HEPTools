import itertools
from copy import deepcopy


import pylhe
import vector 
import awkward as ak


class LHEparse:

    """
    Class for parsing LHE files to lists of 4vectors, awkward arrays and
    eventually ROOT files
    Particles are identified by the pdgid number as defined in the class
    dictionary PDGID
    Args:
    - 
    """
      
    PDGID = {
      1.0:"d"               ,   -1.0:"dbar",
      2.0:"u"               ,   -2.0:"ubar",
      3.0:"s"               ,   -3.0:"sbar",
      4.0:"c"               ,   -4.0:"cbar",
      5.0:"b"               ,   -5.0:"ubar",
      6.0:"t"               ,   -6.0:"tbar",
      7.0:"bprime"          ,   -7.0:"bprimebar",
      8.0:"tprime"          ,   -8.0:"tprimebar",
      11.0:"eP"             ,   -11.0:"eM",
      12.0:"v_eP"           ,   -12.0:"v_eM",
      13.0:"muP"            ,   -13.0:"muM",
      14.0:"v_muP"          ,   -14.0:"v_muM",
      15.0:"tauP"           ,   -15.0:"tauM",
      16.0:"v_tauP"         ,   -16.0:"v_tauM",
      17.0:"tauprimeP"      ,   -17.0:"tauprimeM",
      18.0:"v_tauprimeP"    ,   -18.0:"v_tauprimeM",
      21.0:"g"              , 
      22.0:"A"              , 
      23.0:"Z"              , 
      24.0:"WP"             , "-24":"WM",
      25.0:"h"      
      }
    

    def __init__(self,file_name):
        self.file_name      = file_name
        self.arr            =  pylhe.to_awkward(pylhe.read_lhe_with_attributes("events.lhe"))
        self.list_of_events = []

    def generate_list_of_4vectors(self):

        """
        
        """

        processed_event_template = {v:[] for v in self.PDGID.values()}

        list_of_events = []
        for event in self.arr:
            processed_event = deepcopy(processed_event_template)
            for particle in event["particles"]:
                    processed_event[self.PDGID[particle["id"]]].append(particle["vector"])
            pe = {k:v for k,v in processed_event.items() if len(v)!=0}
            list_of_events.append(pe)

        return list_of_events
    

    def generate_awkward_array(self):
        if len(self.list_of_events) == 0:
            return ak.Array(self.generate_list_of_4vectors())
        else:
            return ak.Array(self.list_of_events)



def parse_LHE_file(filename):
     ParseObject = LHEparse(filename)
     ParseObject.generate_list_of_4vectors()
     return ParseObject.list_of_events 
