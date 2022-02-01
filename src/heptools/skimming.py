# Python3.6 ROOT macro for skimming TTrees


from ROOT import TTree, TFile, vector, gROOT, TLorentzVector


def ResetBranches(dictionary):

	for key in dictionary:
		dictionary[key].resize(0)


def BookBranch(name, suffix, tree):

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
    return branch


def generic_tree_skim(input_tree,branches2keep,**kwargs):

    # Define new tree name
    new_tree_name   = kwargs["tree_name"] if "tree_name" in kwargs else input_tree.GetName() + "_skimmed"
    output_filename = kwargs["file_name"] if "file_name" in kwargs else new_tree_name+".root"
    
    output_file = TFile(output_filename, "RECREATE")

    # Output TTree
    output_tree   = TTree(new_tree_name,   new_tree_name)

    # Turn off all irrelevant branches, and turn back on the pertinent ones
    input_tree.SetBranchStatus("*",0)

    for branch in branches2keep:
        input_tree.SetBranchStatus(branch, 1)

    # Store branches here
    branch_dict = {} 
    for branch in branches2keep:
        branch_dict[branch] = BookBranch(branch, "F", output_tree)        

    total_events = input_tree.GetEntries()

    # Max events
    iteration_cut_off = kwargs["cut_off"] if "cut_off" in kwargs else None

    for counter, event in enumerate(input_tree):

        if counter % 1000 == 0: 
            print("Event",counter,"/",total_events)

        if counter==iteration_cut_off:
            break
            
        ResetBranches(branch_dict) #Don't tell me what to do!
        
        input_tree.GetEntry(counter) #this
        
        for branch in branches2keep:

            if hasattr(event,branch):
                event_branch = getattr(event,branch)

                # If the branch entry is a float or an int
                if type(event_branch)==int or type(event_branch)==float:
                    branch_dict[branch].push_back(event_branch)

                # If the branch entry is a cpp.stl vector
                elif len(event_branch)!=0:
                    event_vec = event_branch
                    for i in event_vec:
                        branch_dict[branch].push_back(i)

            
        output_tree.Fill()

    output_tree.Write()




        

