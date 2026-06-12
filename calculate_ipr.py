# This module calculates the IPR of the localized states of the static disordered Hamiltonian #

import numpy as np

def main(ichain,n_site,folder_path):

    try:
        with open(str(folder_path+'localized_state_coefficient_'+str(ichain)+'.out'),'r') as file:
            filelines=file.readlines()
    except ValueError:
        index=1
        return(index)

    outputfile=open(str(folder_path+'localized_state_ipr_'+str(ichain)+'.out'),'w+')

    for ii in range(n_site):
        c=np.array([float(kk) for kk in filelines[ii].strip().split()])
        ipr=1.0/np.sum(np.power(c,4))
        outputfile.write('{:.10}'.format(ipr)+'\n')

    outputfile.close()

    index=0; return(index)

