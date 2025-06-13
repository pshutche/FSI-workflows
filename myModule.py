import os
import numpy as np
from decimal import *
from tabulate import tabulate

def i(Step, PreviousTimestepSize): # Step is actually a real number!
    
    Skip_steps = 2
    initial_timestep = 1e-2 # this is needed as SyC has nowhere else to get it from and will set 0.0 s otherwise!
    max_timestep = 1e-2 # that which gives 1000 timestep optimum
    target_cfl = 1.0 # 1.0 to 3.0 suggestwed
    max_increase_factor = 2.0 # 1.1 to 2.0 suggested
    Extract_Step = int(Step) - Skip_steps # Fluent report will always be one step behind SyC
    
    if os.path.exists("Fluent_1/max-cfl-rfile.out") == False:
        print("fluent report file doesn't exist, timestep is unchanged..")
        return initial_timestep
    
    elif Extract_Step < Skip_steps:
        print("timestep is unchanged in this step")
        return initial_timestep
    
    else:    
        cflfile = np.genfromtxt("Fluent_1/max-cfl-rfile.out", skip_header = 3, skip_footer = 0) # column 0 = step, 1 = time, 2 = moving mesh cfl, 3 = vof cfl * IMPORTANT: set this order in report def/file in Fluen
        meshcfl = cflfile[Extract_Step,2]
        vofcfl = cflfile[Extract_Step,3]
        #cfl = max(meshcfl,vofcfl)
        cfl = meshcfl
        getcontext().prec = 2 # set number of significant figures
        new_dt = Decimal(PreviousTimestepSize) * min(Decimal(target_cfl) / Decimal(cfl), Decimal(max_increase_factor)) # truncate the timestep
        new_timestep = min(float(new_dt),max_timestep) # convert timestep back to float, limit to max timestep
        c = [["Target CFL",target_cfl],["moving mesh cfl",cfl],["VOF cfl",cflfile[Extract_Step,3]],["max increase factor",max_increase_factor],["max timestep",max_timestep],["old timestep",PreviousTimestepSize],["new timestep",new_timestep]]
        print(tabulate(c))
        return new_timestep