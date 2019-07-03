import time
T = time.time()
from fipy import Variable, FaceVariable, CellVariable, Grid1D, DiffusionTerm,UpwindConvectionTerm,ConvectionTerm, TransientTerm, DiffusionTerm, Viewer,ExponentialConvectionTerm
from fipy.tools import numerix
import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import sys, os
srcPath = os.path.join(os.path.dirname(__file__),'..')
sys.path.append('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/libs')
sys.path.append('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/src/@PhreeqcRM')
sys.path.append('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/src/@IRM_RESULT')
from PhreeqcRM import *

def createPhreeqcInstance(filename,kinetics):

    nx = 100
    ny = 1
    nz = 1
    dx = 1.5
    mesh = Grid1D(nx=nx,dx=dx)
    v_cell = mesh.cellVolumes

    nxyz=nx*ny*nz # number of cells
    nthreads=1 # number of threads
    phreeqc_rm = PhreeqcRM(nxyz, nthreads)
    # create a phreeqc instance
      # Set properties
    status = phreeqc_rm.RM_SetErrorHandlerMode(1)
    status = phreeqc_rm.RM_SetComponentH2O(0)
    status = phreeqc_rm.RM_SetRebalanceFraction(0.5)
    status = phreeqc_rm.RM_SetRebalanceByCell(1)
    status = phreeqc_rm.RM_UseSolutionDensityVolume(0)
    status = phreeqc_rm.RM_SetPartitionUZSolids(0)
    status = phreeqc_rm.RM_SetFilePrefix("Advect_cpp")
    status = phreeqc_rm.RM_OpenFiles()


    status = phreeqc_rm.RM_SetUnitsSolution(2)      # 1, mg/L 2, mol/L 3, kg/kgs
    status = phreeqc_rm.RM_SetUnitsPPassemblage(1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
    status = phreeqc_rm.RM_SetUnitsExchange(1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
    status = phreeqc_rm.RM_SetUnitsSurface(1)       # 0, mol/L cell 1, mol/L water 2 mol/L rock
    status = phreeqc_rm.RM_SetUnitsGasPhase(1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
    status = phreeqc_rm.RM_SetUnitsSSassemblage(1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
    status = phreeqc_rm.RM_SetUnitsKinetics(1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
      # Set conversion from seconds to days
      # I prefer to work with seconds, so I keep it to 1.0
    status = phreeqc_rm.RM_SetTimeConversion(1.0/86400)
    # Set representative volume
    rv = np.ones_like(v_cell,dtype='float64') # volume is in liter
    status = phreeqc_rm.RM_SetRepresentativeVolume(rv)

      # Set initial porosity
    por = np.ones(nxyz,dtype='float64')*1
    status = phreeqc_rm.RM_SetPorosity(por)
     # Set initial saturation
    sat = np.ones(nxyz,dtype='float64')
    status = phreeqc_rm.RM_SetSaturation(sat)
    	# Set cells to print chemistry when print chemistry is turned on
    print_chemistry_mask = np.zeros(nxyz,dtype='int')
    print_chemistry_mask[:np.int(nxyz/2)]=1
    status = phreeqc_rm.RM_SetPrintChemistryMask(print_chemistry_mask)
    	# Partitioning of uz solids
    #status = phreeqc_rm.RM_SetPartitionUZSolids(0)
    	# Demonstation of mapping, no symmetry, only one row of cells
    grid2chem = np.arange(nxyz,dtype='int32')

    status = phreeqc_rm.RM_CreateMapping(grid2chem)
    if (status < 0):
        status = phreeqc_rm.RM_DecodeError(status)
    nchem = phreeqc_rm.RM_GetChemistryCellCount()
    status = phreeqc_rm.RM_SetPrintChemistryOn(0, 1, 0) # workers, initial_phreeqc, utility
    	# Set printing of chemistry file
    status = phreeqc_rm.RM_LoadDatabase("minteq.v4.dat")

    	# Demonstration of error handling if ErrorHandlerMode is 0
      # commented out for now (AAE)
    	# if (status != Iphreeqc_rm.RM_OK)
    	# {
    	# 	l = phreeqc_rm.RM_GetErrorStringLength()
    	# 	errstr = (char *) malloc((size_t) (l * sizeof(char) + 1))
    	# 	phreeqc_rm.RM_GetErrorString(errstr, l+1)
    	# 	fprintf(stderr,"Beginning of error string:\n")
    	# 	fprintf(stderr,"%s", errstr)
    	# 	fprintf(stderr,"End of error string.\n")
    	# 	free(errstr)
    	# 	errstr = NULL
    	# 	phreeqc_rm.RM_Destroy()
    	# 	exit(1)
    	# }
    	# Run file to define solutions and reactants for initial conditions, selected output
    	# There are three types of IPhreeqc instances in Phreeqcphreeqc_rm.RM
    	# Argument 1 refers to the workers for doing reaction calculations for transport
    	# Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
    	# Argument 3 refers to the Utility instance
    status = phreeqc_rm.RM_RunFile(1, 1, 1, filename)
    	# Clear contents of workers and utility
    Str="DELETE -all"
    status = phreeqc_rm.RM_RunString(1, 0, 1, Str)	# workers, initial_phreeqc, utility
    	# Determine number of components to transport
    ncomps = phreeqc_rm.RM_FindComponents()

    gfw = np.zeros(ncomps,dtype='float64')
    status = phreeqc_rm.RM_GetGfw(gfw)

    #for i in range(ncomps):
    # I consider the maximum number of characters for components' name to be 10
    components = np.zeros(ncomps,dtype='U10')
    for i in range(ncomps):
        status = phreeqc_rm.RM_GetComponent(i, components,10)
    #    print(components[i].encode(), "\t", gfw[i], "\n")


    ic1 = -1*np.ones(nxyz*7,dtype=np.int32)
    ic2 = -1*np.ones(nxyz*7,dtype=np.int32)
    f1 = np.ones(nxyz*7,dtype=np.float64)
    for i in np.arange(nxyz):
        ic1[i] = 1       # Solution 1
        ic1[i+nxyz] = -1      # Equilibrium phases none
        ic1[i+2*nxyz] = -1       # Exchpdange 1
        ic1[i+3*nxyz]= -1      # Surface none
        ic1[i+4*nxyz]= -1      # Gas phase none
        ic1[i+5*nxyz]= -1      # Solid solutions none
        ic1[i+6*nxyz]= kinetics      # Kinetics none
    status = phreeqc_rm.RM_InitialPhreeqc2Module(ic1, ic2, f1);
    t = 0.0 # changed time to time1 to avoid conflicts with julia time function
    time_step = 0.0
    c = np.zeros(ncomps*nxyz,dtype=np.float64)
    status = phreeqc_rm.RM_SetTime(t)
    status = phreeqc_rm.RM_SetTimeStep(time_step)
    status = phreeqc_rm.RM_SetSpeciesSaveOn(1)

    status = phreeqc_rm.RM_RunCells()
    status = phreeqc_rm.RM_GetConcentrations(c)

    nbound = 1
    bc1 = np.zeros(nbound,dtype='int32')
    bc2 = -1* np.ones(nbound,dtype='int32')
    bc_f1 = np.ones(nbound,dtype='float64')
    bc_conc = np.zeros(ncomps*nbound,dtype='float64')
    status = phreeqc_rm.RM_InitialPhreeqc2Concentrations(bc_conc, nbound, bc1, bc2, bc_f1)
    density = np.ones(nxyz,dtype=np.float64)
    volume = np.zeros(nxyz,dtype=np.float64)
    pressure = 2*np.ones(nxyz,dtype=np.float64)
    temperature = 20*np.ones(nxyz,dtype=np.float64)
    sat_calc = np.zeros(nxyz,dtype=np.float64)

    status = phreeqc_rm.RM_SetDensity(density)
    status = phreeqc_rm.RM_SetPressure(pressure)
    status = phreeqc_rm.RM_SetTemperature(temperature)






    return bc_conc,mesh,ncomps,components,c,nxyz,phreeqc_rm
