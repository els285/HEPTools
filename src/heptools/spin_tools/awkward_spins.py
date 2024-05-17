import awkward as ak
import pylhe
import vector 
import numpy as np
import sys
import time 

def parse(filename: str):

    """
    Parses the LHE and generates an awkward array of the top, anti-top and
    charged leptons Momentum4D vectors
    """
    array = pylhe.to_awkward(pylhe.read_lhe_with_attributes('events.lhe'))

    # Pull out the relevant particle 4-vectors
    top            = array.particles[array.particles.id == 6] ["vector"]
    anti_top       = array.particles[array.particles.id ==-6] ["vector"]
    pos_elec       = array.particles[array.particles.id ==-11]["vector"]
    neg_elec       = array.particles[array.particles.id ==11] ["vector"]
    pos_mu         = array.particles[array.particles.id ==-13]["vector"]
    neg_mu         = array.particles[array.particles.id ==13] ["vector"]

    leptonP = ak.concatenate([pos_elec,pos_mu],axis=1)
    leptonM = ak.concatenate([neg_elec,neg_mu],axis=1)

    out = ak.zip({
            "eventinfo": array.eventinfo,
            "tops": top,
            "anti_tops": anti_top,
            "pos_leptons": leptonP,
            "neg_leptons": leptonM
        }, depth_limit=1, with_name="Event")
    
    return  out 


def boost(arr: "pylhe.awkward.EventArray"):

    """
    Performs the relevant boosts and scalar products to return the cos-variables
    """

    top      = arr["tops"]
    anti_top = arr["anti_tops"]
    leptonP  = arr["pos_leptons"]
    leptonM  = arr["neg_leptons"]

    # Build ttbar system
    ttbar = top + anti_top 

    # Boost tops lab --> ttbar_CoM
    top_in_CoM      = top.boostCM_of(ttbar)
    antitop_in_CoM  = anti_top.boostCM_of(ttbar) 

    # Boost leptons lab --> ttbar_CoM --> parent_tops'_CoM
    leptonP_in_CoM     = leptonP.boostCM_of(ttbar)
    leptonM_in_CoM     = leptonM.boostCM_of(ttbar)
    leptonP_in_top     = leptonP_in_CoM.boostCM_of(top_in_CoM)
    leptonM_in_antitop = leptonM_in_CoM.boostCM_of(antitop_in_CoM)

    # Compute the unit 3-vector directions
    Kdirection        = top_in_CoM.to_beta3().unit()
    leptonP_direction = leptonP_in_top.to_beta3().unit()
    leptonM_direction = leptonM_in_antitop.to_beta3().unit()

    return Kdirection, leptonP_direction, leptonM_direction



def helicity_basis_observables(Kdirection, leptonP_direction, leptonM_direction):

    """
    The helicity basis is a basis used in the measurement of spin observables at
    the LHC
    """

    # Kinematics
    z     = vector.obj(x=0,y=0,z=1)
    cos_T = z.dot(Kdirection)
    sin_T = (1 - cos_T**2)**0.5
    mask = 1*(cos_T > 0) -1*(cos_T < 0)

    # Helicity basis observables
    cos_K_plus  = Kdirection.dot(leptonP_direction)
    cos_K_minus = -Kdirection.dot(leptonM_direction) 

    Ndirection = z.cross(Kdirection)/sin_T
    cos_N_plus  = mask*Ndirection.dot(leptonP_direction)
    cos_N_minus = -mask*Ndirection.dot(leptonM_direction) 

    Rdirection  =  1/sin_T * (z - cos_T*Kdirection)
    cos_R_plus  = mask*Rdirection.dot(leptonP_direction)
    cos_R_minus = -mask*Rdirection.dot(leptonM_direction) 

    cos_phi     = leptonP_direction.dot(leptonM_direction)

    helicity_observables = ak.zip({
    "cos_K_plus"  :  cos_K_plus,
    "cos_K_minus" :  cos_K_minus,
    "cos_N_plus"  :  cos_N_plus,
    "cos_N_minus" :  cos_N_minus,
    "cos_R_plus"  :  cos_R_plus,
    "cos_R_minus" :  cos_R_minus,
    "cos_phi"     :  cos_phi 
    }, depth_limit=1, with_name="Event")

    return helicity_observables


def compute_spin_parameters(observable_array):

    """
    Computes the spin parameters (C_{i,j}, B_k^{+/-})
    """

    Ckk = -9*(np.multiply(observable_array["cos_K_plus"][:,0].to_numpy() , observable_array["cos_K_minus"][:,0].to_numpy()).mean())
    Cnn = -9*(np.multiply(observable_array["cos_N_plus"][:,0].to_numpy() , observable_array["cos_N_minus"][:,0].to_numpy()).mean())
    Crr = -9*(np.multiply(observable_array["cos_R_plus"][:,0].to_numpy() , observable_array["cos_R_minus"][:,0].to_numpy()).mean())
    Crk = -9*(np.multiply(observable_array["cos_R_plus"][:,0].to_numpy() , observable_array["cos_K_minus"][:,0].to_numpy()).mean())
    Ckr = -9*(np.multiply(observable_array["cos_K_plus"][:,0].to_numpy() , observable_array["cos_R_minus"][:,0].to_numpy()).mean())
    Cnr = -9*(np.multiply(observable_array["cos_N_plus"][:,0].to_numpy() , observable_array["cos_R_minus"][:,0].to_numpy()).mean())
    Crn = -9*(np.multiply(observable_array["cos_R_plus"][:,0].to_numpy() , observable_array["cos_N_minus"][:,0].to_numpy()).mean())
    Cnk = -9*(np.multiply(observable_array["cos_N_plus"][:,0].to_numpy() , observable_array["cos_K_minus"][:,0].to_numpy()).mean())
    Ckn = -9*(np.multiply(observable_array["cos_K_plus"][:,0].to_numpy() , observable_array["cos_N_minus"][:,0].to_numpy()).mean())

    CrkP = Crk + Ckr 
    CrkM = Crk - Ckr 
    CnrP = Cnr + Crn
    CnrM = Cnr - Crn 
    CnkP = Cnk + Ckn 
    CknM = Cnk - Ckn

    BkP = -3*observable_array["cos_K_plus"][:,0].to_numpy().mean()
    BkM = -3*observable_array["cos_K_minus"][:,0].to_numpy().mean()
    BnP = -3*observable_array["cos_N_plus"][:,0].to_numpy().mean()
    BnM = -3*observable_array["cos_N_minus"][:,0].to_numpy().mean()
    BrP = -3*observable_array["cos_R_plus"][:,0].to_numpy().mean()
    BrM = -3*observable_array["cos_R_minus"][:,0].to_numpy().mean()

    return {"Ckk"       : Ckk,
            "Cnn"       : Cnn,
            "Crr"       : Crr,
            "CrkP"      : CrkP,
            "CrkM"      : CrkM,
            "CnrP"      : CnrP,
            "CnrM"      : CnrM,
            "CnkP"      : CnkP,
            "CknM"      : CknM,
            "BkP"       : BkP,
            "BkM"       : BkM,
            "BnP"       : BnP,
            "BnM"       : BnM,
            "BrP"       : BrP,
            "BrM"       : BrM}


def histograms(observable_array, Nbins:int=10 ):

    """
    Use the cosvariable arrays to build the relevant angular observable
    histograms
    Args:
    - Observable_array: awkward-array of the cos-variables
    - Nbins: integer defining the binning (optional)
    """

    import boost_histogram as bh

    DH = {  "Ckk"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "Cnn"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "Crr"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "CrkP"      : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "CrkM"      : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "CnrP"      : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "CnrM"      : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "CnkP"      : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "CknM"      : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "BkP"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "BkM"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "BnP"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "BnM"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "BrP"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "BrM"       : bh.Histogram(bh.axis.Regular(Nbins, -1, +1)),
            "cos_phi"   : bh.Histogram(bh.axis.Regular(Nbins, -1, +1))
    }

    # Diagonal terms
    DH["Ckk"].fill(np.multiply(observable_array["cos_K_plus"][:,0].to_numpy() , observable_array["cos_K_minus"][:,0].to_numpy()))
    DH["Cnn"].fill(np.multiply(observable_array["cos_N_plus"][:,0].to_numpy() , observable_array["cos_N_minus"][:,0].to_numpy()))
    DH["Crr"].fill(np.multiply(observable_array["cos_R_plus"][:,0].to_numpy() , observable_array["cos_R_minus"][:,0].to_numpy()))

    # Cross-terms
    Crk = np.multiply(observable_array["cos_R_plus"][:,0].to_numpy() , observable_array["cos_K_minus"][:,0].to_numpy())
    Ckr = np.multiply(observable_array["cos_K_plus"][:,0].to_numpy() , observable_array["cos_R_minus"][:,0].to_numpy())
    Cnr = np.multiply(observable_array["cos_N_plus"][:,0].to_numpy() , observable_array["cos_R_minus"][:,0].to_numpy())
    Crn = np.multiply(observable_array["cos_R_plus"][:,0].to_numpy() , observable_array["cos_N_minus"][:,0].to_numpy())
    Cnk = np.multiply(observable_array["cos_N_plus"][:,0].to_numpy() , observable_array["cos_K_minus"][:,0].to_numpy())
    Ckn = np.multiply(observable_array["cos_K_plus"][:,0].to_numpy() , observable_array["cos_N_minus"][:,0].to_numpy())

    DH["CrkP"].fill(np.add(Crk,Ckr))
    DH["CrkM"].fill(np.add(Crk,-Ckr))
    DH["CnrP"].fill(np.add(Cnr,Crn))
    DH["CnrM"].fill(np.add(Cnr,-Crn))
    DH["CnkP"].fill(np.add(Cnk,Ckn))
    DH["CknM"].fill(np.add(Cnk,-Ckn))

    # Polarisations
    DH["BkP"].fill(observable_array["cos_K_plus"][:,0].to_numpy())
    DH["BkM"].fill(observable_array["cos_K_minus"][:,0].to_numpy())
    DH["BnP"].fill(observable_array["cos_N_plus"][:,0].to_numpy())
    DH["BnM"].fill(observable_array["cos_N_minus"][:,0].to_numpy())
    DH["BrP"].fill(observable_array["cos_R_plus"][:,0].to_numpy())
    DH["BrM"].fill(observable_array["cos_R_minus"][:,0].to_numpy())

    DH["cos_phi"].fill(observable_array["cos_phi"][:,0].to_numpy())

    return DH

