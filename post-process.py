import os
import numpy as np
import math
from numpy.lib.stride_tricks import sliding_window_view
import matplotlib.pyplot as plt
from tabulate import tabulate

# Collection of post-processing tools for FSI analysis
# Can parse file.nlh from MAPDL, and report definitions from Fluent
# Time dependent files
# Generate FFT
# Extract max frequency from FFT
# Extract max displacement from file.nlh
# Create time and frequency plots
# Create frequency vs velocity plot (multi-run, or velocity sweep for flutter etc...)
# Create displacement vs velocity plot (multi-run, or velocity sweep for flutter etc...)
# Ansys is not responsible for results obtained using this python script
# This script is delivered as an example of what can be done in python
# Script could be combined with pySystemCoupling/pyFluent/pyMAPDL scripts or called from them

# Global parameters - can be put into .json file eventually
Analyse_Run = 3
Total_Runs = 10
First_velocity = 2.0 # m/s
Velocity_increment = 0.5 # m/s
Last_velocity = First_velocity + Velocity_increment * (Total_Runs - 1) # m/s
Run_velocity = First_velocity + (Analyse_Run - 1) * 0.5 # m/s
f_n = 9.1 # Hz
Omega_n = 2 * math.pi * f_n
D = 3 * 0.0254 # inches -> m
Domain_length = 83 * 0.0254 # inches -> m
FTT = Domain_length / Run_velocity # s
Freq_res = 0.5 # Hz
Timestep = 0.001 # s
Sampling_frequency = 1 / Timestep
Nyquist_max = Sampling_frequency / 2
Num_samples = 1 / Freq_res / Timestep
Sampling_time = Num_samples * Timestep
Start_time_sampling = 5 * FTT # s
End_time_sampling = Start_time_sampling + Sampling_time # s

# Time window for averaging
time_window = Sampling_time #s
window_size = int(time_window / Timestep) # last elements to consider at a time
#window_size = 5000

clfile = np.genfromtxt("Run"+str(Analyse_Run)+"/cl-rfile.out", skip_header = 3, skip_footer = 0)
cdfile = np.genfromtxt("Run"+str(Analyse_Run)+"/cd-rfile.out", skip_header = 3, skip_footer = 0)
cl = np.array(clfile[-window_size:,1]) # keep 2nd column (variable) and consider sampling time window
cd = np.array(cdfile[-window_size:,1]) # keep 2nd column (variable) and consider sampling time window
cd = np.array(cdfile[-window_size:,1]) # keep 2nd column (variable)
time = clfile[-window_size:,2] # keep 1st column (time) and consider sampling time window

# Extract displacements from file.nlh, removing duplicate time points
dispfile = np.genfromtxt("Run"+str(Analyse_Run)+"/file.nlh", skip_header = 7, skip_footer = 2) # needs cleaning of repeated column 0 rows
# First separate into 2 columns
first, second = np.array(dispfile[:,0]), np.array(dispfile[:,1])
# Unique rows and indices from first array
unique = np.unique(first, return_index = True, return_counts = True)
# Obtain unique times, indices of the first of unique times and unique counts all as separate arrays
uniques, first_indices, counts = unique
# Find indices of last uniques
last_indices = first_indices + (counts-1)
# Extract these indices from original array (1D array is OK)
dispfull = np.array(second[last_indices,])
# Apply time window
disp = np.array(dispfull[-window_size:])

# Time averaged CD
cdavg = np.average(cd)

# Max displacement (mm)
maxdisp = max(abs(disp)) * 1e3

# Calculate FFT for vortex shedding frequency

n = cl.size
fourier = np.fft.fft(cl)
fourierabs = np.abs(fourier)/n
freq = np.fft.fftfreq(n, d = Timestep)
fvs = round(abs(freq[np.argmax(fourierabs),]),2)

# Calculate FFT for frequency of flow induced vibration

n = disp.size
fourier = np.fft.fft(disp)
fourierabs = np.abs(fourier)/n
freq = np.fft.fftfreq(n, d = Timestep)
viv = round(abs(freq[np.argmax(fourierabs),]),2)

# Print run characteristics
c = [["Run no.",Analyse_Run], ["Inlet velocity (m/s)",Run_velocity], ["Domain length (m)",Domain_length], ["Flow through time (s)",FTT], ["Freq resolution (Hz)",Freq_res], ["Timestep (s)",Timestep], ["Number of samples",Num_samples], 
["Sampling start time (s)",Start_time_sampling], ["Sampling end time (s)",End_time_sampling], ["Time-averaged CD",cdavg], ["Vortex shedding frequency (Hz)",fvs], ["VIV frequency (Hz)",viv], ["Max displacement (mm)",maxdisp]]
print(tabulate(c))

fig, [ax1, ax2, ax3, ax4] = plt.subplots(nrows = 4, ncols = 1)
ax1.plot(time,cl, '.-')
ax2.plot(time,cd, '.-')
ax3.plot(time,disp, '.-')
ax4.plot(freq,fourierabs, '.-')

ax1.set_xlabel("Time (s)")
ax2.set_xlabel("Time (s)")
ax3.set_xlabel("Time (s)")
ax4.set_xlabel("Frequency (Hz)")

ax1.set_ylabel("Lift coefficient")
ax2.set_ylabel("Drag coefficient")
ax3.set_ylabel("Displacement (m)")

ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

#ax1.set_xlim(0,int(End_time_sampling))
#ax2.set_xlim(0,int(End_time_sampling))
#ax3.set_xlim(0,int(End_time_sampling))
#ax4.set_xlim(0,Nyquist_max)
ax4.set_xlim(0,20)

#ax1.set_ylim(-1,1)
#ax1.set_ylim(0.2,1.2)
#ax3.set_ylim(0.85,0.95)
#ax4.set_ylim(0.85,0.95)

plt.tight_layout()
plt.show()

# Populate table with input velocities
dispmatrix = np.zeros((Total_Runs,2))
#freqmatrix = np.zeros((Total_Runs,2))
freqmatrix = np.zeros((Total_Runs,3)) 
velocities = np.arange(First_velocity,Last_velocity+1e-8,0.5)
run = np.arange(1,Total_Runs+1)
#dispmatrix[:,0] = velocities
dispmatrix[:,0] = velocities/Omega_n/D
#freqmatrix[:,0] = velocities
freqmatrix[:,0] = velocities/Omega_n/D # fill in column 0 (1st column)

i = 0
while i <= Total_Runs-1:
    run[i] = i + 1
# Extract displacements from file.nlh, removing duplicate time points
    if os.path.exists("Run"+str(run[i])+"/file.nlh") == False:
        print("file Run"+str(run[i])+" not found, proceeding to next run...")
        dispmatrix[i,1] = None
        freqmatrix[i,1] = None
        run[i] = 0
        i = i + 1
        continue # go to next loop
    else:
        print("file Run"+str(run[i])+" found")
        dispfile = np.genfromtxt("Run"+str(run[i])+"/file.nlh", skip_header = 7, skip_footer = 2) # needs cleaning of repeated column 0 rows
# First separate into 2 columns
        first, second = np.array(dispfile[:,0]), np.array(dispfile[:,1])
# Unique rows and indices from first array
        unique = np.unique(first, return_index = True, return_counts = True)
# Obtain unique times, indices of the first of unique times and unique counts all as separate arrays
        uniques, first_indices, counts = unique
# Find indices of last uniques
        last_indices = first_indices + (counts-1)
# Extract these indices from original array (1D array is OK)
        dispfull = np.array(second[last_indices,])
# Apply time window
        disp = np.array(dispfull[-window_size:])
        maxdisp = max(abs(disp))
        #dispmatrix[i,1] = maxdisp # 2nd column becomes displacement
        dispmatrix[i,1] = maxdisp/D # 2nd column becomes normalised displacement
# Cl
        clfile = np.genfromtxt("Run"+str(run[i])+"/cl-rfile.out", skip_header = 3, skip_footer = 0)
        cl = np.array(clfile[-window_size:,1]) # keep 2nd column (variable) and consider sampling time window
# Calculate FFT of cl
        n = cl.size
        fourier = np.fft.fft(cl) # change input signal here
        fourierabs = np.abs(fourier)/n
        freq = np.fft.fftfreq(n, d = Timestep)
        vs = round(abs(freq[np.argmax(fourierabs),]),2)
        #freqmatrix[i,1] = vs # 2nd column becomes vs frequency
        freqmatrix[i,1] = vs/f_n # 2nd column becomes vs frequency
# Calculate FFT of disp
        n = disp.size
        fourier = np.fft.fft(disp) # change input signal here
        fourierabs = np.abs(fourier)/n
        freq = np.fft.fftfreq(n, d = Timestep)
        viv = round(abs(freq[np.argmax(fourierabs),]),2)
        #freqmatrix[i,2] = viv # 3rd column becomes viv frequency
        freqmatrix[i,2] = viv/f_n # 3rd column becomes viv frequency
        i = i + 1

# Velocity sweep!
    
min_v = min(velocities)
max_v = max(velocities)
min_f = 0.198 * min_v / D / f_n
max_f = 0.198 * max_v / D / f_n
st_line = np.array([[1,min_f],[Total_Runs,max_f]]) # ([fmin,umin], [fmax,umax])
print(st_line)

fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1)

#ax1.plot(run,freqmatrix[:,1], '.-',run,freqmatrix[:,2], '.-')
ax1.plot(run,freqmatrix[:,1],'o',run,freqmatrix[:,2],'o',st_line[:,0],st_line[:,1],'-')
ax2.plot(run,dispmatrix[:,1],'o')

ax1.grid()
ax2.grid()

ax1.set_xlabel("Run")
ax2.set_xlabel("Run")

ax1.set_ylabel("f/f_n")
ax2.set_ylabel("y/D")

ax1.set_xticks(run)
ax2.set_xticks(run)

ax1.legend(['freq_vs', 'freq_st', 'St = 0.198'], loc = 'upper left')

def velfromrun(run):
    return (First_velocity + (run - 1) * 0.5) / Omega_n / D
    
def runfromvel(run):
    return (First_velocity + (run - 1) * 0.5)  / Omega_n / D

secax = ax2.secondary_xaxis("top", functions = (velfromrun, runfromvel))
secax.set_xlabel("U = v / (Omega_n * D)")

plt.tight_layout()
plt.show()