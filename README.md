The project contains code written in Python and Fortran 90.

**User Environment:**
To successfully run the code, the user must have the following activated in the user environment
-- python=3.12 
-- numpy=2.4.4

**Compilation:**
The Fortran 90 codes have been compiled using numpy.f2py tool using Intel Fortran compiler. The user should import the 
Intel compiler module before the compilation. The compilation code is given in the file **f2py_compile**.

**Running the code:**
The primary code the user runs is **polymer_intra_chain_mobility.py**. The user can choose to use the **DEFAULT_INPUT** 
Or make a new input maintaining the format. If no argument is provided for the input filename, the code will use the default 
input. 

**python polymer_intra_chain_mobility.py** 
Or
**python polymer_intra_chain_mobility.py 'Input File Name'**

**Input:** 
The format for the input file is as follows: 

**1. comment line  
2. blank line  
3. 3. n_chain: No. of instances of the randomized 1D model Hamiltonian of the polymer; must be an integer and >=1  
4. property: Physical property to be determined; must be a string - either 'mobility' or 'diffusivity' or 'ipr' or 'll'  
5. n_site: No. of SRU in polymer chain; must be an integer and >=2  
6. obc_or_pbc: Open/periodic boundary condition applied to the polymer chain; must be a string - either 'obc' or 'pbc'  
7. alpha: Average on-site energy  
8. beta: Average electronic coupling  
9. sigma_alpha_static: Standard deviation in the static disorder of on-site energy  
10. sigma_beta_static: Standard deviation in static disorder of the coupling  
11. rate_equation_type: Type of rate equation expression to be used; must be a string - either 'marcus' or 'jortner' or 'm-a' or 'generalized'  
12. sigma_alpha_dynamic: Dynamic disorder in on-site energy  
13. sigma_beta_dynamic: Dynamic disorder in coupling  
14. lambda: Reorganization energy/SRU due to coupling with vibrational modes  
15. T: Temperature  
16. field: Electric field strength  
17. sru_length: SRU length  
18. carrier: carrier type; must be a string - either 'h' or 'e'  
19. blank line**  

**Output:**
The program currently calculates the energies of localized states of the disordered Hamiltonian, their coefficients, their 
inverse participation ratio, their centroid position, localization length, hopping rates between them and their steady-state
occupation probability. The output files are stored in the **Output_Folder** folder.

The mobility list is written in the **mobility_list.dat** file in the **main (runtime) folder**.
