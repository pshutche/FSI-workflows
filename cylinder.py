import os
import ansys.fluent.core as pyfluent
import ansys.mapdl.core as pymapdl
import ansys.systemcoupling.core as pysystemcoupling
from tabulate import tabulate
import json
import shutil
import numpy as np

Fluent_folder = "Fluent_1"      
#Fluent_case = "Fluent.cas.h5"
Fluent_case = "Fluent-overset-laminar.cas.h5"
Fluent_cores = 32
#Fluent_ui_mode = "no_gui_or_graphics" # -g
Fluent_ui_mode = "no_gui" # -gu
#Fluent_data = "Fluent.RANS.h5"
Mechanical_folder = "Mechanical_2"
Mechanical_jobname = "Mechanical"
Mechanical_input = "Mechanical.dat"
Mechanical_scp = "Mechanical_2/Mechanical.scp"
Mechanical_cores = 16
Domain_length = 83 * 0.0254 # m
Freq_res = 0.5 # Hz
Timestep = 0.001 # s
Num_samples = 1 / Freq_res / Timestep # = number of timesteps
Total_Runs = 10
run = np.arange(1,Total_Runs+1)
First_velocity = 2.0 # m/s
Velocity_increment = 0.25 # m/s
Last_velocity = First_velocity + Velocity_increment * (Total_Runs - 1) # m/s
velocities = np.arange(First_velocity, Last_velocity + 1e-8, Velocity_increment)

# Start Mechanical
#mechanical = pymapdl.launch_mapdl(nproc = Mechanical_cores,run_location = Mechanical_folder,override = True,jobname = Mechanical_jobname)
mechanical = pymapdl.launch_mapdl(nproc = Mechanical_cores,override = True)

# Stream input file (not working, bug logged)
#mechanical.input(fname = Mechanical_input)
#mechanical.etlist() # test the model is loaded

# Start Fluent
fluent = pyfluent.launch_fluent(case_file_name = Fluent_case, processor_count = Fluent_cores, ui_mode = Fluent_ui_mode, cwd = Fluent_folder, graphics_driver = "null")

# pyfluent datamodel shortcuts
fluset = fluent.settings

# Fluent batch options
fluidbo = fluset.file.batch_options
fluidbo.confirm_overwrite = False 
# (default = True)
fluidbo.exit_on_error = True 
# (default = False)
fluidbo.hide_answer = True
# (default = False)

# MDM controls
sp = fluset.setup.dynamic_mesh.methods.smoothing
#sp.method = "radial"
#sprad = sp.radial_settings
#sprad.smooth_bl_with_adj()
#sprad.smooth_bl_with_adj = True
#sp.method = "diffusion"
#spdif = sp.diffusion_settings
#spdif.diffusion_coeff_function = "boundary-distance"
#spdif.smooth_from_ref = "True"
#spdif.diffusion_coeff_parameter = 1.5
#spdif.boundary_distance_method = "False"
#spdif.max_iter = 100
#spdif.relative_tolerance = 1e-7
#spdif.verbosity = 2

#sp.method = "spring"
#spspr = sp.spring_settings
#spspr.max_iter = 500
#spspr.constant_factor = 0.1
#spspr.tolerance = 1e-4
#spspr.verbosity = 2

# Initialise Fluent
Inlet_velocity = velocities[0]
fluset.setup.reference_values.velocity = Inlet_velocity
fluidin = fluset.solution.initialization
fluidin.initialization_type = "standard"
fluidin.defaults = {'x-velocity': Inlet_velocity}
fluidin.initialize()

# Start System Coupling
syc = pysystemcoupling.launch(start_output=True)

# Connect participants
#solid = syc.setup.add_participant(participant_session = mechanical)
solid = syc.setup.add_participant(input_file = Mechanical_scp)
fluid = syc.setup.add_participant(participant_session = fluent)

# Create interface
fsi_interface = syc.setup.add_interface(
    side_one_participant=solid,side_one_regions=['FSIN_1'],
    side_two_participant=fluid,side_two_regions=['fsi'])

# Create data transfers
force_transfer_name = syc.setup.add_data_transfer(
    interface=fsi_interface,
    target_side="One",
    target_variable="FORC",
    source_variable="force",
)

displacement_transfer_name = syc.setup.add_data_transfer(
    interface=fsi_interface,
    target_side="Two",
    source_variable="INCD",
    target_variable="displacement",
)

# Apply transformation if necessary

# Set timestep
syc.setup.solution_control.time_step_size = Timestep

# Start parametric sweep

i = 0
while i <= Total_Runs - 1:
    run[i] = i + 1
    Inlet_velocity = velocities[i]
    FTT = Domain_length / Inlet_velocity # s
    Start_time_sampling = 5 * FTT # s
    End_time_sampling = Start_time_sampling + (Num_samples * Timestep) # s
    if i == 0:
        end_time = End_time_sampling
    else:
        end_time = end_time + End_time_sampling
# Print run characteristics
    c = [["Inlet velocity (m/s)",Inlet_velocity], ["Domain length (m)",Domain_length], ["Flow through time (s)",FTT], ["Freq resolution (Hz)",Freq_res], ["Timestep (s)",Timestep], ["Number of samples",Num_samples], ["Sampling start time (s)",Start_time_sampling], ["Sampling end time (s)",End_time_sampling]]
    print(tabulate(c))

# Set inlet velocity (uncorrected)
    fluset.setup.boundary_conditions.velocity_inlet["inlet"].momentum.velocity.value = Inlet_velocity

# Solve
    #steps = int(End_time_sampling/Timestep)    
    #syc.setup.solution_control.duration_option = "NumberOfSteps"
    #syc.setup.solution_control.duration_option = "EndTime"
    syc.setup.solution_control.end_time = end_time
    #syc.setup.solution_control.number_of_steps = 100000 # arbitrarily large
    #syc.solution.step(count = steps)
    syc.solution.solve()

# Move reports and file.nlh to its permanent home
    directory_name = "Run"+str(run[i])
    shutil.move(Mechanical_folder + "/" + "file.nlh", directory_name + "/" + "file.nlh")
    shutil.move(Fluent_folder + "/" + "cd-rfile.out", directory_name + "/" + "cd-rfile.out")
    shutil.move(Fluent_folder + "/" + "cl-rfile.out", directory_name + "/" + "cl-rfile.out")

# End parametric sweep

syc.solution.shutdown()
fluent.exit()
#mechanical.exit()
#syc.exit()
