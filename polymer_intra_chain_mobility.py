# MAIN CODE FOR THE CALCULATION OF CHARGE MOBILITY OR CHARGE DIFFUSIVITY IN BULK POLYMER FROM A GENERALIZED 1D MODEL HAMILTONIAN 

# The user can choose -- 
#	A. charge mobility within standard masters equation approach ('mobility_steady_state')
# 	B. charge mobility from diffusivity via master equation + Einstein relation ('mobility_diffusivity')

# The Input file -- 
#	I. default input: DEFAULT_INPUT
#  II. user-defined input: filename must be provided by the user

# There are 16 paramaters in the input file which are organized in terms of their use in the modules (see main text). 
# Energy, length and electric field units are electronvolt (eV), angstrom and volts/cm respectively.

# 	Main program:
#	1. n_chain: No. of instances of the randomized 1D model Hamiltonian of the polymer; must be an integer and >=1
#	2. property: Physical property to be determined; must be a string - either 'mobility_steady_state' or 'mobility_diffusivity' (optionally with 'ipr', 'll')


#	Module 1:
#	3. n_site: No. of SRU in polymer chain (for the definition of SRU, please see the references in the main text); must be an integer and >=200 for steady_state_mobility, >600 for diffusivity_mobility.
#	4. obc_or_pbc: Open/periodic boundary condition applied to the polymer chain; must be a string - either 'obc' or 'pbc'
#	5. alpha: Average on-site energy
#	6. beta: Average electronic coupling
#	7. sigma_alpha_static: Standard deviation in the static disorder of on-site energy
#	8. sigma_beta_static: Standard deviation in static disorder of the coupling 


#	Module 2 and 3:
#	9. rate_equation_type: Type of rate equation expression to be used (please see the main document for references); must be a string - either 'marcus' or 'jortner' or 'm-a' or 'generalized'
#   10. sigma_alpha_dynamic: Dynamic disorder in on-site energy 
#   11. sigma_beta_dynamic: Dynamic disorder in coupling
#   12. lambda: Reorganization energy/SRU due to coupling with vibrational modes
#   13. T: Temperature
#   14. field: Electric field strength
#   15. sru_length: SRU length
#   16. carrier: carrier type; must be a string - either 'h' or 'e'

# INPUT FILE FORMAT SPECFICATION: 
#   1ST LINE: STARTS WITH A # SYMBOL
#   2ND LINE: KEEP IT BLANK
#   LAST LINE: KEEP IT BLANK
#   KEEP A BLANK SPACE BETWEEN MAIN INPUT AND OPTIONAL INPUT SECTION


# ******************************************************************************************************************************************************************************************** #

import sys
import random
import numpy as np
import os
from pathlib import Path
import check_input
import static_disordered_hamiltonian
import calculate_ipr
import calculate_ll
import hopping_rate
import steady_state_mobility 
import diffusivity
import diffusivity3d

def main():

# Define Inputfile

    if len(sys.argv) == 1 :
        print('Default inputfile is used')
        inputfilename='DEFAULT_INPUT'
    elif len(sys.argv) == 2 :
        print('User-defined inputfile is used')
        inputfilename=sys.argv[1]
    else:
        print('Error in command line execution - too many arguments; exit calculation')
        sys.exit(1)
 

# Check Inputfile for formatting error

    check_flag=check_input.main(inputfilename)
    if check_flag != 0 :
        print('Error in the format of input file; exit calculation')
        sys.exit(1)


# Read Inputfile

    with open(inputfilename,'r') as inputfile:
        lines=[line.strip() for line in inputfile]
    if len(lines) <= 21:
        n_chain=int(lines[2])
        property=lines[3].strip().split() 
        n_site=int(lines[4])
        obc_or_pbc=lines[5]
        alpha=float(lines[6])
        beta=float(lines[7])
        sigma_alpha_static=float(lines[8])
        sigma_beta_static=float(lines[9])
        rate_equation_type=lines[10]
        sigma_alpha_dynamic=float(lines[11])
        sigma_beta_dynamic=float(lines[12])
        Lambda=float(lines[13])
        T=float(lines[14])
        field=float(lines[15])
        sru_length=float(lines[16])
        carrier=lines[17]

        print('Main input reading - done')

    output_folder_path=Path('Output_Folder')
    output_folder_path.mkdir(exist_ok=True)
    output_folder = str(output_folder_path.resolve()) + '/'

 
# Define the fixed seed for random number generation and for reproducibility purposes; can be modified by user if required.

    random_seed=n_chain+37*n_site+41*int(T)
    random.seed(random_seed)

    if any(string in property for string in ('mobility_steady_state','mobility_diffusivity', 'mobility_3d' )):
        mobility=np.zeros((n_chain))

# Loop over the number of randomized polymer conformations

    for ichain in range(1,n_chain+1):

        print('Random Hamiltonian -',ichain,flush=True)
        random_int=random.randint(1,int(1.0e9))

#       Call module 1 for the construction of static disordered Hamiltonian

        hamiltonian_index=static_disordered_hamiltonian.hamiltonian(np.int32(ichain),np.int32(n_site),
                        np.float64(alpha),np.float64(beta),np.float64(sigma_alpha_static),
                        np.float64(sigma_beta_static),np.int32(random_int),obc_or_pbc,output_folder)

###        print(static_disordered_hamiltonian.hamiltonian.__doc__)

        if hamiltonian_index != 0:
            print('Error - Hamiltonian construction and/or diagonlaization',flush=True); sys.exit(1)
        else:
            print('Localized states are calculated for random Hamiltonian',flush=True)

        
        if 'ipr' in property:

            ipr_index=calculate_ipr.main(ichain,n_site,output_folder)        

            if ipr_index != 0:
                print('Error - IPR calculation',flush=True)
            else:
                print('IPR of localized states are calculated',flush=True)


        if any(string in property for string in ('ll','mobility_steady_state','mobility_diffusivity','mobility_3d')):

            ll_index=calculate_ll.main(ichain,n_site,sru_length,obc_or_pbc,output_folder)        

#            print(calculate_ll.main.__doc__)

            if ll_index != 0:
                print('Error - localization length calculation',flush=True)
            else:
                print('Localization length of localized states are calculated',flush=True)


        if any(string in property for string in ('mobility_steady_state','mobility_diffusivity','mobility_3d')) and ll_index == 0:

            if rate_equation_type == 'marcus':

                rate_index=hopping_rate.marcus(ichain,n_site,sru_length,sigma_alpha_dynamic,
                            sigma_beta_dynamic,Lambda,T,field,carrier,obc_or_pbc,output_folder)

#                print(hopping_rate.marcus.__doc__)

                if rate_index != 0:
                    print('Error - hopping rate calculation',flush=True)
                else:
                    print('Hopping rate between localized states are calculated',flush=True)


        if 'mobility_steady_state' in property and rate_index == 0:

                mobility[ichain-1],mobility_index=steady_state_mobility.main(ichain,n_site,sru_length,field,carrier,obc_or_pbc,output_folder)

#                print(steady_state_mobility.main.__doc__)

                if mobility_index != 0:
                    print('Error - mobility calculation',flush=True)
                else:
                    print('Steady state intra-chain mobility is calculated',flush=True)


        if 'mobility_diffusivity' in property and rate_index == 0:

                mobility[ichain-1],mobility_index=diffusivity.main(ichain,n_site,sru_length,carrier,T,output_folder)

#                print(diffusivity.main.__doc__)

                if mobility_index != 0:
                    print('Error - mobility calculation',flush=True)
                else:
                    print('Intra-chain mobility is calculated from diffusivity',flush=True)

    if 'mobility_3d' in property:
        V_inter = float(os.environ.get('V_INTER', '0.03'))
        P = float(os.environ.get('P_PERSIST', '50.0'))

        D_3D, mu_3D, index_3d = diffusivity3d.main( n_chain, n_site, sru_length, carrier, T, output_folder, V_inter=V_inter, P=P) 

        if index_3d != 0:
            print('Error - 3D diffusivity calculation', flush=True)
        else:
            print(f'3D mobility = {mu_3D:.6e} cm^2/(V.s)', flush=True)
    
    if any(string in property for string in ('mobility_steady_state','mobility_diffusivity','mobility_3d' )):
        with open('mobility_list.dat','w+') as mobility_file:
            for ii in range(0,n_chain):
                mobility_file.write(str('{:.10f}'.format(mobility[ii]))+'\n')
        print('Calculation of mobility for all random Hamiltonian is complete',flush=True)
 

main()
