#!/usr/bin/env python
'''
Updated by Ethan Simpson January 28th 2021.
Changes:
    - Writing to TTree during parsing, rather than looping over events twice. Large reduction in running time
    - Split into separate functions, including the "main" function, called "parser", which can be called from another python macro.
'''

# Python modules
import argparse
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from ROOT import TTree, TFile, vector, gROOT, TLorentzVector, TVector3
from copy import copy
import math
import time

from math import log, tan, acos, pi, copysign

###--- Argument parser setup ---###
parser = argparse.ArgumentParser()
parser.add_argument("--infile",  
                    help    = "File in LHE format to be read",  
                    nargs   = '?', 
                    type    = str, 
                    default = "random.lhe")
parser.add_argument("--outfile", 
                    help    = "File in ROOT format for output", 
                    nargs   = '?', 
                    type    = str, 
                    default = "test.root")
parser.add_argument("--outtree", 
                    help    = "output TTree name",              
                    nargs   = '?', 
                    type    = str, 
                    default = "test")
# parser.add_argument("--cme",
#                     help    = "Centre of mass energy in MEV, default is 13 TeV",
#                     nargs   = '?',
#                     type    = float,
#                     default = 13000
#                         )   
args = parser.parse_args()

include_decays = True


class particle(object):
    def __init__(self, index, pid, status, mother1, mother2, color1, color2, px, py, pz, e, mass, f1, spin):
        self.index    = index
        self.pid      = pid
        self.status   = status
        self.mother1  = mother1
        self.mother2  = mother2
        self.color1   = color1
        self.color2   = color2
        self.px       = px
        self.py       = py
        self.pz       = pz
        self.e        = e
        self.m        = mass
        self.f1       = f1
        self.spin     = spin

        # Computed on the fly
        self.p2      = self.px**2 + self.py**2 + self.pz**2
        self.p       = self.p2 ** 0.5
        self.pt      = (self.px**2 + self.py**2)**0.5
        self.theta   = acos(self.pz / self.p)
        try:
            self.phi = copysign( acos(self.px / self.pt), self.py )
        except ZeroDivisionError:
            self.phi = 0.
        try:
            self.eta = -log(tan(self.theta / 2.))
        except ValueError:
            self.eta = copysign(9999.0, self.pz / self.p)
        try:
            self.y   = 0.5*log((self.e + self.pz)/(self.e - self.pz))
        except (ZeroDivisionError, ValueError):
            self.y   = 0

        self.helicity = 0 #3 = left, -3 = right
        if self.spin * self.pz > 0:
            self.helicity = -3
        else:
            self.helicity = 3

def deltaR(a, b):
    deta = a.eta - b.eta
    dphi = a.phi - b.phi
    if dphi > pi:
        dphi -= pi
    if dphi < -pi:
        dphi += pi
    return (deta**2. + dphi**2.)**.5

def ResetBranches(branches):

    for branch in branches.values():
        branch.resize(0)


def BookBranch(name, suffix, tree, branches):

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
    branches[name]=branch
    return branch


def cos_helicity(top, parent_top, ttbar, lepton, sign):

    """ Function to calculate the helicity angle
    top        = the top qaurk 4 vector
    parent_top = the parent t or tbar of the lepton 4 vector
    ttbar      = the 4 vector of the ttbar system (i.e. t + tbar)
    lepton     = the 4 vector fo the lepton
    sign       = a necessary sign argument (comes from the theory, basically
                 it's related to the charge of the lepton)

    WARNING: the boost functions here will change the 4 vectors you're passing
    to this function. Always pass copies if you intend to use them again later
    in the code.
    """

    ###-- First, we need to define a boost vector, because we're going to
    ###-- boost our lepton into, first the ttbar frame, and then the top
    ###-- rest frame.
    boost_to_ttbar = ttbar.BoostVector()
    boost_to_ttbar = boost_to_ttbar*-1.

    ###-- Boost everything into the ttbar system
    parent_top.Boost(boost_to_ttbar)
    top.Boost(boost_to_ttbar)
    lepton.Boost(boost_to_ttbar)

    ###-- Now work out how to boost to the top rest frame
    boost_to_parent_top = parent_top.BoostVector()
    boost_to_parent_top = boost_to_parent_top*-1.


    ###-- Boost the lepton there
    lepton.Boost(boost_to_parent_top)

    ###-- Now work out the spin analysing basis (which were calling the k vec)
    k_vector = top.Vect().Unit()
    k_vector = k_vector*sign

    ###-- Now get the angle between the lepton and kvector and return it
    cos_theta = lepton.Vect().Unit().Dot(k_vector)


    if (math.isnan(cos_theta)):
        return -55.
    else:
        return cos_theta

def cos_transverse(top, parent_t, ttbar, lepton, sign):

    boost_to_ttbar = ttbar.BoostVector()
    boost_to_ttbar *= -1.
    
    parent_t.Boost(boost_to_ttbar)
    top.Boost(boost_to_ttbar)
    lepton.Boost(boost_to_ttbar)
    
    boost_to_parent_t = parent_t.BoostVector()
    boost_to_parent_t *= -1.

    lepton.Boost(boost_to_parent_t)
  
    k_vector = top.Vect().Unit()
    #k_vector *= sign
    p_vector = TVector3(0,0,1)

    y = p_vector.Dot(k_vector)
    r = pow((1. - y*y),0.5)

    n_vector = (1./r)*(p_vector.Cross(k_vector));

    if (sign == 1):
        if(y > 0):
            n_vector *= 1.
        if(y < 0):
            n_vector *= -1.
    elif (sign == -1):
        if(y > 0):
            n_vector *= -1.
        if(y < 0):
            n_vector *= 1.

    theta = lepton.Vect().Unit().Dot(n_vector)

    if (math.isnan(theta)):
        return -55.
    else:
        return theta


def cos_raxis(top, parent_t, ttbar, lep, sign):

    boost_to_ttbar = ttbar.BoostVector()
    boost_to_ttbar *= -1.
      
    parent_t.Boost(boost_to_ttbar)
    top.Boost(boost_to_ttbar)
    lep.Boost(boost_to_ttbar)
      
    boost_to_parent_t = parent_t.BoostVector()
    boost_to_parent_t *= -1.    
    lep.Boost(boost_to_parent_t)
    
    k_vector = top.Vect().Unit()
    #k_vector *= sign
    p_vector = TVector3(0,0,1)  
    y = p_vector.Dot(k_vector)
    r = pow((1. - y*y),0.5) 
    r_vector = (1./r)*(p_vector - y*k_vector)

    if (sign == 1):
        if(y > 0):
            r_vector *= 1.
        if(y < 0):
            r_vector *= -1.
    elif (sign == -1):
        if(y > 0):
            r_vector *= -1.;
        if(y < 0):
            r_vector *= 1.;

    theta = lep.Vect().Unit().Dot(r_vector)

    if (math.isnan(theta)):
        return -55.
    else:
        return theta



def initialise_branches(outtree):

    '''
    Specified list of branches to be made on ROOT output TTree
    Initialises each branch which is a ROOT vector
    Generates the branches dictionary stores all these branches
    '''
    branches = dict()

    init_branches = [
    'init_pz',
    'init_pdgid',
    'init_hel',
    'init_spin'
    ]
    
    weight_branches = ['weights','reweight1']

    production_branches = [
    'top_pt',
    'top_eta',
    'top_phi',
    'top_e',
    'top_m',
    'top_y',
    'tbar_pt',
    'tbar_eta',
    'tbar_phi',
    'tbar_e',
    'tbar_m',
    'tbar_y',
    'ttbar_pt',
    'ttbar_eta',
    'ttbar_phi',
    'ttbar_e',
    'ttbar_m',
    'ttbar_y',
    'G_pt',
    'G_eta',
    'G_phi',
    'G_e',
    'G_m',
    'G_spin',
    'onshell_top',
    'onshell_tbar']

    decay_branches = [
    'cosp_hel',
    'cosm_hel',
    'cosp_raxis',
    'cosm_raxis',
    'cosp_trans',
    'cosm_trans',
    'b_pt',
    'b_eta',
    'b_phi',
    'b_e',
    'b_m',
    'b_spin',
    'bbar_pt',
    'bbar_eta',
    'bbar_phi',
    'bbar_e',
    'bbar_m',
    'bbar_spin',
    'Wp_pt',
    'Wp_eta',
    'Wp_phi',
    'Wp_e',
    'Wp_m',
    'Wm_pt',
    'Wm_eta',
    'Wm_phi',
    'Wm_e',
    'Wm_m',
    'lp_pt',
    'lp_eta',
    'lp_phi',
    'lp_e',
    'lp_m',
    'lp_pdgid',
    'lp_spin',
    'lm_pt',
    'lm_eta',
    'lm_phi',
    'lm_e',
    'lm_m',
    'lm_pdgid',
    'lm_spin',
    'dphi_ll',
    'v_pt',
    'v_eta',
    'v_phi',
    'v_e',
    'v_m',
    'v_pdgid',
    'v_spin',
    'vbar_pt',
    'vbar_eta',
    'vbar_phi',
    'vbar_e',
    'vbar_m',
    'vbar_pdgid',
    'vbar_spin']


    for branch_name in init_branches:
        BookBranch(branch_name,'/F', outtree, branches)

    for branch_name in weight_branches:
        BookBranch(branch_name,'/F', outtree, branches)

    for branch_name in production_branches:
        if branch_name == 'onshell_top' or branch_name == 'onshell_tbar':
            BookBranch(branch_name,'/I', outtree, branches)
        else:
            BookBranch(branch_name,'/F', outtree, branches)

    for branch_name in decay_branches:
        BookBranch(branch_name,'/F', outtree, branches)

    return branches

def push2ROOT(outtree,particlelist,branches):

    '''
    LHE event particle ---> ROOT output
        --input: particlelist = a list of particles parsed from LHE file
        Pushes all particle attributes to outtree branches
    '''

    # Particle assignment
    hasTop=False
    hasTbar=False
    has_lp   = False
    has_lm   = False


    for counter, p in enumerate(particlelist):
        if counter < 2:
            branches["init_pz"].push_back(p.pz)
            branches["init_pdgid"].push_back(p.pid)
            branches["init_hel"].push_back(p.helicity) 
            branches["init_spin"].push_back(p.spin)               

        if p.pid == 6:
            branches["top_pt"].push_back(p.pt)
            branches["top_eta"].push_back(p.eta)
            branches["top_phi"].push_back(p.phi)
            branches["top_e"].push_back(p.e)
            branches["top_m"].push_back(p.m)
            branches["top_y"].push_back(p.y)
            hasTop = True
        elif p.pid == -6:
            branches["tbar_pt"].push_back(p.pt)
            branches["tbar_eta"].push_back(p.eta)
            branches["tbar_phi"].push_back(p.phi)
            branches["tbar_e"].push_back(p.e)
            branches["tbar_m"].push_back(p.m)
            branches["tbar_y"].push_back(p.y)
            hasTbar = True
        elif p.pid == 5:
            branches["b_pt"].push_back(p.pt)
            branches["b_eta"].push_back(p.eta)
            branches["b_phi"].push_back(p.phi)
            branches["b_e"].push_back(p.e)
            branches["b_m"].push_back(p.m)
            branches["b_spin"].push_back(p.spin)                               
        elif p.pid == -5:
            branches["bbar_pt"].push_back(p.pt)
            branches["bbar_eta"].push_back(p.eta)
            branches["bbar_phi"].push_back(p.phi)
            branches["bbar_e"].push_back(p.e)
            branches["bbar_m"].push_back(p.m)
            branches["bbar_spin"].push_back(p.spin)                                               
        elif p.pid == 24:
            branches["Wp_pt"].push_back(p.pt)
            branches["Wp_eta"].push_back(p.eta)
            branches["Wp_phi"].push_back(p.phi)
            branches["Wp_e"].push_back(p.e)
            branches["Wp_m"].push_back(p.m)
        elif p.pid == -24:
            branches["Wm_pt"].push_back(p.pt)
            branches["Wm_eta"].push_back(p.eta)
            branches["Wm_phi"].push_back(p.phi)
            branches["Wm_e"].push_back(p.e)
            branches["Wm_m"].push_back(p.m)
        elif p.pid == 11 or p.pid == 13 or p.pid == 15:
            branches["lm_pt"].push_back(p.pt)
            branches["lm_eta"].push_back(p.eta)
            branches["lm_phi"].push_back(p.phi)
            branches["lm_e"].push_back(p.e)
            branches["lm_m"].push_back(p.m)
            branches["lm_pdgid"].push_back(p.pid)
            branches["lm_spin"].push_back(p.spin)   
            has_lm = True                                                            
        elif p.pid == -11 or p.pid == -13 or p.pid == -15:
            branches["lp_pt"].push_back(p.pt)
            branches["lp_eta"].push_back(p.eta)
            branches["lp_phi"].push_back(p.phi)
            branches["lp_e"].push_back(p.e)
            branches["lp_m"].push_back(p.m)
            branches["lp_pdgid"].push_back(p.pid)
            branches["lp_spin"].push_back(p.spin)  
            has_lp = True                                                             
        elif p.pid == 12 or p.pid == 14 or p.pid == 16:
            branches["v_pt"].push_back(p.pt)
            branches["v_eta"].push_back(p.eta)
            branches["v_phi"].push_back(p.phi)
            branches["v_e"].push_back(p.e)
            branches["v_m"].push_back(p.m)
            branches["v_pdgid"].push_back(p.pid)
            branches["v_spin"].push_back(p.spin)                                                               
        elif p.pid == -12 or p.pid == -14 or p.pid == -16:
            branches["vbar_pt"].push_back(p.pt)
            branches["vbar_eta"].push_back(p.eta)
            branches["vbar_phi"].push_back(p.phi)
            branches["vbar_e"].push_back(p.e)
            branches["vbar_m"].push_back(p.m)
            branches["vbar_pdgid"].push_back(p.pid)
            branches["vbar_spin"].push_back(p.spin)                                                               
        elif abs(p.pid) == 21:
            branches["G_pt"].push_back(p.pt)
            branches["G_eta"].push_back(p.eta)
            branches["G_phi"].push_back(p.phi)
            branches["G_e"].push_back(p.e)
            branches["G_m"].push_back(p.m)
            branches["G_spin"].push_back(p.spin)                                                               
        #else:
            #print("LHEparser: WARNING - Particle",p.pid,"wasn't assigned to anything!")

    ###-- Do ttbar --###
    top   = TLorentzVector()
    tbar  = TLorentzVector()
    if hasTop:
        top.SetPtEtaPhiM(branches["top_pt"][0],  branches["top_eta"][0],  branches["top_phi"][0],  branches["top_m"][0])
        branches["onshell_top"].push_back(1)
    else:
        _lepton   = TLorentzVector()
        _neutrino = TLorentzVector()
        _bottom   = TLorentzVector()
        _lepton.SetPtEtaPhiM(   branches["lp_pt"][0], branches["lp_eta"][0], branches["lp_phi"][0], branches["lp_m"][0])
        _neutrino.SetPtEtaPhiM( branches["v_pt"][0],  branches["v_eta"][0],  branches["v_phi"][0],  branches["v_m"][0])
        _bottom.SetPtEtaPhiM(   branches["b_pt"][0],  branches["b_eta"][0],  branches["b_phi"][0],  branches["b_m"][0])
        top = _lepton + _neutrino + _bottom
        branches["onshell_top"].push_back(0)
        branches["top_pt"].push_back(top.Pt())
        branches["top_eta"].push_back(top.Eta())
        branches["top_phi"].push_back(top.Phi())
        branches["top_e"].push_back(top.E())
        branches["top_m"].push_back(top.M())
        branches["top_y"].push_back(top.Y())

    if hasTbar:
        tbar.SetPtEtaPhiM(branches["tbar_pt"][0], branches["tbar_eta"][0], branches["tbar_phi"][0], branches["tbar_m"][0])
        branches["onshell_tbar"].push_back(1)
    else:
        _lepton   = TLorentzVector()
        _neutrino = TLorentzVector()
        _bottom   = TLorentzVector()
        _lepton.SetPtEtaPhiM(   branches["lm_pt"][0],   branches["lm_eta"][0],   branches["lm_phi"][0],   branches["lm_m"][0])
        _neutrino.SetPtEtaPhiM( branches["vbar_pt"][0], branches["vbar_eta"][0], branches["vbar_phi"][0], branches["vbar_m"][0])
        _bottom.SetPtEtaPhiM(   branches["bbar_pt"][0], branches["bbar_eta"][0], branches["bbar_phi"][0], branches["bbar_m"][0])
        tbar = _lepton + _neutrino + _bottom
        branches["onshell_tbar"].push_back(0)
        branches["tbar_pt"].push_back(tbar.Pt())
        branches["tbar_eta"].push_back(tbar.Eta())
        branches["tbar_phi"].push_back(tbar.Phi())
        branches["tbar_e"].push_back(tbar.E())
        branches["tbar_m"].push_back(tbar.M())
        branches["tbar_y"].push_back(tbar.Y())
    ttbar = top + tbar

    branches["ttbar_pt"].push_back(ttbar.Pt())
    branches["ttbar_eta"].push_back(ttbar.Eta())
    branches["ttbar_phi"].push_back(ttbar.Phi())
    branches["ttbar_e"].push_back(ttbar.E())
    branches["ttbar_m"].push_back(ttbar.M())
    branches["ttbar_y"].push_back(ttbar.Rapidity())

    if include_decays:

        if has_lp and has_lm:
            lep_p = TLorentzVector()
            lep_m = TLorentzVector()        

            lep_p.SetPtEtaPhiM(branches["lp_pt"][0],  branches["lp_eta"][0],  branches["lp_phi"][0],  branches["lp_m"][0])
            lep_m.SetPtEtaPhiM(branches["lm_pt"][0],  branches["lm_eta"][0],  branches["lm_phi"][0],  branches["lm_m"][0])   

            branches["dphi_ll"].push_back(lep_p.DeltaPhi(lep_m))

            cos_hel_p   = cos_helicity(  copy(top), copy(top),  copy(ttbar), copy(lep_p), +1)
            cos_hel_m   = cos_helicity(  copy(top), copy(tbar), copy(ttbar), copy(lep_m), -1)        
            cos_trans_p = cos_transverse(copy(top), copy(top),  copy(ttbar), copy(lep_p), +1)
            cos_trans_m = cos_transverse(copy(top), copy(tbar), copy(ttbar), copy(lep_m), -1) 
            cos_raxis_p = cos_raxis(     copy(top), copy(top),  copy(ttbar), copy(lep_p), +1)
            cos_raxis_m = cos_raxis(     copy(top), copy(tbar), copy(ttbar), copy(lep_m), -1) 
        
            branches["cosp_hel"].push_back(cos_hel_p)
            branches["cosm_hel"].push_back(cos_hel_m)
            branches["cosp_raxis"].push_back(cos_raxis_p)
            branches["cosm_raxis"].push_back(cos_raxis_m)
            branches["cosp_trans"].push_back(cos_trans_p)
            branches["cosm_trans"].push_back(cos_trans_m)

    # outtree.Fill()
    
def parse_event(elem,branches,outtree):

    '''
    The actual parsing...
        --input: an element of the XML file which corresponds to an interaction events
        Strips the strings and constructs event_data plus particle objects

    '''

    LHEparticlelist = elem.text.split("\n")
    LHEparticlelist = filter(None, LHEparticlelist)
    LHEparticlelist = list(LHEparticlelist)

    event_data = [float(value) for value in LHEparticlelist[0].split()]

    num_particles   = event_data[0]
    pid             = event_data[1]
    nominal_weight  = event_data[2]
    fac_scale       = event_data[3]
    aqed            = event_data[4]
    aqcd            = event_data[5]

    nominal_weight = float(LHEparticlelist[0].split()[2])
    branches["weights"].push_back(nominal_weight)

    LHEparticlelist.pop(0)

    particlelist = []

    '''
    Looping over rest of data in event
    Generating particles
    Pushing to ROOT TTree branches
    '''

    for iLHEp, LHEp in enumerate(LHEparticlelist):

        ###--- cleanup the particles ---###
        if "pdf" in LHEp: #Skip lines with pdf info, need to store this. Not needed for LHEF3
            continue 
        if "#aMCatNLO" in LHEp: #Skip weird lines
            continue

        LHEp = LHEp.strip()
        LHEp = LHEp.split()
        if LHEp == []:
            continue

        try: 
            LHEp = list(map(int, LHEp[:6])) + list(map(float, LHEp[6:]))
        except:
            print("Failed to convert list: ",LHEp)
            raise AssertionError

        assert len(LHEp)==13, "Failed Assertion"


        particlelist.append( particle(iLHEp, *LHEp) )

    ### Call the push2ROOT function which will write each particle's parameters to the outtree
    push2ROOT(outtree,particlelist,branches)

    # Fill the TTree with this event's particles and weight
    outtree.Fill()


def parser(infile,outfile,outtree):


    print("LHEparser: Welcome to the LHE Parser, attempting to parse input file")
    print("LHEparser: (this may take a while if you have a lot of events)")   

    t0 = time.time()

    #Initialise output ROOT file
    gROOT.cd()
    outfile = TFile(outfile, "recreate")
    outtree = TTree(outtree, outtree)
    gROOT.cd()

    # Initialise branches
    branches=initialise_branches(outtree)

    ### Parse the XML 
    assert not ".lhe.gz" in infile, "You need to un-tar the input file first"

    counter =0 
    for (index, elem) in ET.iterparse(infile):

        ResetBranches(branches)
        if elem.tag == "event":
            counter += 1
            if counter % 1000 == 0: print("LHEparser: Event number ",str(counter))
            parse_event(elem,branches,outtree)

        elem.clear()

    outfile.Write()
    outfile.Close()

    print(time.time()-t0)


def main():

    parser(infile=args.infile,outfile=args.outfile,outtree=args.outtree)

if __name__ == '__main__':
    main()

    # exit()


