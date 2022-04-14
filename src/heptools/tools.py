def compute_ratio(ROOT_hist_numer,ROOT_hist_denom):

    d_hist = ROOT_hist_numer.Clone()
    d_hist.Divide(ROOT_hist_numer,ROOT_hist_denom)

    return d_hist


def print_key_names(file):
    [print(k.GetName()) for k in file.GetListOfKeys()]

def print_branch_names(ttree):
    [print(k.GetName()) for k in ttree.GetListOfBranches()]    