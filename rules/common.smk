#Create contrasts lists
def get_contrasts(conditions):
    #Remove underscores from conditions to facilitate DeSeq2 contrasts in R
    #Doesn't have to use underscores (actually may not want to)
    underscore_flag = any("_" in cond for cond in conditions) 
    if underscore_flag:
        print("WARNING: Some 'Conditions' in sample manifest contains underscores '_' : Changing to hyphens '-'")
        conditions = [cond.replace("_", "-") for cond in conditions]

    contrasts = list(itertools.combinations(set(conditions), 2))
    return(sorted(['_'.join(map(str,sorted(pair))) for pair in contrasts]))
