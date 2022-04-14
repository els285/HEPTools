def key_names(file):
    [print(k.GetName()) for k in file.GetListOfKeys()]

def branch_names(ttree):
    [print(k.GetName()) for k in ttree.GetListOfBranches()]
