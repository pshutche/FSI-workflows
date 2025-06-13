# FSI-workflows
Useful python scripts for FSI workflows
by Paul Hutcheson, at some point in 2025 when he had some free time :)

-----------
cylinder.py
-----------
Does a parametric sweep changing Fluent inlet velocity
Outputs cd, cl and MAPDL deflection against time to separate runN folders for post-processing (that can be done consequently using post-process.py)

-----------
myModule.py
-----------
Adaptive time-stepping module to be hooked into SyC GUI expression functions
Note: there is an active bug in versions pre 2026R1 that results in changes to the module code not updating in the data model. Instead a previous version is used.
Workaround: It is necessary to rename at least the function inside the module. A more robust workaround is to delete the expression function in SyC GUI, . If function name does not appear, reload expression function.
Time-step is CFL driven. A report definition and report file must be created in Fluent with step, time, mesh-based cfl, vof cfl (4 columns). User can choose to take average or max of these two. If vof cfl is used it is recommended to define the report definition with the following expression:

Maximum(ElementConvectionCourantNumber*4*Volumefraction(phase='liq')*Volumefraction(phase='gas'),['fluid'])

where 'fluid' is the name of the cell zone

-----------
post-process.py
-----------
Loops over a series of runs/run folders
Plot time histories of Fluent CD, CL and structural deflection (from standard MAPDL deformation tracker file.nlh) and analyse one FSI run
Performs FFT to extract frequency spectrum and peak frequency over a frequency range
Extracts max displacement over time
Uses a time window for post-processing (e.g. where flow development is needed before sampling begins)
Plots structural and vortex shedding frequencies against velocity for flutter and stability analysis
