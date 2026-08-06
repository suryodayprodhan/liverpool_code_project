# Check for formatting error in the Input

import numpy as np

def is_integer(s):
    try:
        int(s.strip())
        return True
    except ValueError:
        return False

def is_float(s):
    try:
        float(s.strip())
        return True
    except ValueError:
        return False

def main(filename):

    with open(filename,'r') as inputfile:
        lines=inputfile.readlines()

    if(len(lines) < 19):
        check_flag=1
        return(check_flag)

    integer_index_list=[2,4]
    float_index_list=list(range(6,10))+list(range(11,17))

    for ii in integer_index_list:
        check_int=is_integer(lines[ii])
        if not check_int:
            check_flag=1; return(check_flag)
        if int(lines[ii]) <= 0:
            check_flag=1; return(check_flag)
    
    if any(string not in ('mobility_steady_state','mobility_diffusivity','ll','ipr') for string in lines[3].strip().split()):
        check_flag=1; return(check_flag)
    if lines[5].strip() not in ('obc','pbc'):
        check_flag=1; return(check_flag)
    if lines[10].strip() not in ('marcus','jortner','m-a','generalized'):
        check_flag=1; return(check_flag)
    if lines[17].strip() not in ('h','e'):
        check_flag=1; return(check_flag)

    for ii in float_index_list:
        check_float=is_float(lines[ii])
        if not check_float:
            check_flag=1; return(check_flag)

    for ii in range(2,len(float_index_list)):
        if np.sign(float(lines[float_index_list[ii]])) == -1:
            check_flag=1; return(check_flag)

#    if(len(lines) > 21): Check the same for optional arguments

    check_flag=0
    return(check_flag)

