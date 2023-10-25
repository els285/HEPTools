import awkward as ak
import pylhe
import vector 

def parse_and_boost(filename):

    array = pylhe.to_awkward(pylhe.read_lhe_with_attributes('events.lhe'))

    # Pull out the relevant particle 4-vectors
    top            = array.particles[array.particles.id == 6]["vector"]
    anti_top       = array.particles[array.particles.id==-6]["vector"]
    pos_elec       = array.particles[array.particles.id==11]["vector"]
    neg_elec       = array.particles[array.particles.id==-11]["vector"]
    pos_mu         = array.particles[array.particles.id==13]["vector"]
    neg_mu         = array.particles[array.particles.id==-13]["vector"]

    leptonP = ak.concatenate([pos_elec,pos_mu],axis=1)
    leptonM = ak.concatenate([neg_elec,neg_mu],axis=1)

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
    Kdirection = top_in_CoM.to_beta3().unit()
    leptonP_direction = leptonP_in_top.to_beta3().unit()
    leptonM_direction = leptonM_in_antitop.to_beta3().unit()

    return Kdirection, leptonP_direction, leptonM_direction



def helicity_basis_observables(Kdirection, leptonP_direction, leptonM_direction):

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

    helicity_observables = ak.zip({
    "cos_K_plus"  :  cos_K_plus,
    "cos_K_minus" :  cos_K_minus,
    "cos_N_plus"  :  cos_N_plus,
    "cos_N_minus" :  cos_N_minus,
    "cos_R_plus"  :  cos_R_plus,
    "cos_R_minus" :  cos_R_minus
    }, depth_limit=1, with_name="Event")

    return helicity_observables



