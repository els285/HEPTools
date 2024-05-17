import vector 


""""
Non--ROOT-based tools for computing particular angular observables from which
the spin density matrix elements can be extracted

"""

def boost_to_parent_tops(top:'MomentumObject4D',anti_top:'MomentumObject4D',
             leptonP:'MomentumObject4D',leptonM:'MomentumObject4D'):
    
    """
    Takes MomentumObject4D vectors corresponding to the top-quark,
    anti-top-quark, positively-charged lepton and negatively-charged lepton
    Boosts all 4-vectors into the ttbar CoM frame
    Boosts the leptons in their respective parent top's frame
    Computes the unit 3-vectors of the top in the ttbar CoM, and the leptons in
    their parent tops' frame
    """

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
    direction = top_in_CoM.to_beta3().unit()
    leptonP_direction = leptonP_in_top.to_beta3().unit()
    leptonM_direction = leptonM_in_antitop.to_beta3().unit()

    return direction,leptonP_direction,leptonM_direction


def cos_k(Kdirection,leptonP_direction,leptonM_direction):

    cos_plus  = Kdirection.dot(leptonP_direction)
    cos_minus = -Kdirection.dot(leptonM_direction) 
    return cos_plus, cos_minus


def cos_n(Kdirection,leptonP_direction,leptonM_direction):

    z     = vector.obj(x=0,y=0,z=1)
    cos_T = z.dot(Kdirection)
    sin_T = (1 - cos_T**2)**0.5
    Ndirection = z.cross(Kdirection)/sin_T

    A = 1 if cos_T > 0 else -1

    cos_plus  = A*Ndirection.dot(leptonP_direction)
    cos_minus = -A*Ndirection.dot(leptonM_direction) 
    return cos_plus, cos_minus


def cos_r(Kdirection,leptonP_direction,leptonM_direction):

    z     = vector.obj(x=0,y=0,z=1)
    cos_T = z.dot(Kdirection)
    sin_T = (1 - cos_T**2)**0.5
    Rdirection =  1/sin_T * (z - cos_T*Kdirection)

    A = 1 if cos_T > 0 else -1

    cos_plus  = A*Rdirection.dot(leptonP_direction)
    cos_minus = -A*Rdirection.dot(leptonM_direction) 
    return cos_plus, cos_minus



def helicity_basis_angular_observables(top,antitop,lepP,lepM):

    """
    Builds the helicity-basis angular observables plus cos_phi
    """

    Kdirection, lepP_direction, lepM_direction = boost_to_parent_tops(top,antitop,lepP,lepM)
    obs = {}
    obs["cos_k_p"], obs["cos_k_m"] = cos_k(Kdirection,lepP_direction,lepM_direction)
    obs["cos_n_p"], obs["cos_n_m"] = cos_n(Kdirection,lepP_direction,lepM_direction)
    obs["cos_r_p"], obs["cos_r_m"] = cos_r(Kdirection,lepP_direction,lepM_direction)
    obs["cos_phi"]                 = lepP_direction.dot(lepM_direction)
    return obs


def cos_phi(top,antitop,lepP,lepM):

    Kdirection, lepP_direction, lepM_direction = boost_to_parent_tops(top,antitop,lepP,lepM)
    return lepP_direction.dot(lepM_direction)


def lab_frame_cos_phi(lepP,lepM):
    return lepP.to_beta3().unit().dot(lepM.to_beta3().unit())