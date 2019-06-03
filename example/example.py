import numpy as np
import sys, os
srcPath = os.path.join(os.path.dirname(__file__),'..')
sys.path.append('../src')
from PhreeqcRM import *

Nx=40
Ny=1
Nz=1
nxyz=Nx*Ny*Nz # number of cells
nthreads=1 # number of threads
id=RM_Create(nxyz, nthreads) # create a phreeqc instance
  # Set properties
status = RM_SetErrorHandlerMode(id, 1)
status = RM_SetComponentH2O(id, 0)
status = RM_SetRebalanceFraction(id,0.5)
status = RM_SetRebalanceByCell(id, 1)
status = RM_UseSolutionDensityVolume(id, 0)
status = RM_SetPartitionUZSolids(id, 0)
status = RM_SetFilePrefix(id, "Advect_cpp")
status = RM_OpenFiles(id)


status = RM_SetUnitsSolution(id, 2)      # 1, mg/L 2, mol/L 3, kg/kgs
status = RM_SetUnitsPPassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
status = RM_SetUnitsExchange(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
status = RM_SetUnitsSurface(id, 1)       # 0, mol/L cell 1, mol/L water 2 mol/L rock
status = RM_SetUnitsGasPhase(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
status = RM_SetUnitsSSassemblage(id, 1)  # 0, mol/L cell 1, mol/L water 2 mol/L rock
status = RM_SetUnitsKinetics(id, 1)      # 0, mol/L cell 1, mol/L water 2 mol/L rock
  # Set conversion from seconds to days
  # I prefer to work with seconds, so I keep it to 1.0
status = RM_SetTimeConversion(id, 1.0/86400)
# Set representative volume
rv = np.ones(nxyz,dtype='float64') # volume is in liter
status = RM_SetRepresentativeVolume(id, rv)

  # Set initial porosity
por = np.ones(nxyz,dtype='float64')*0.2
status = RM_SetPorosity(id, por)
 # Set initial saturation
sat = np.ones(nxyz,dtype='float64')
status = RM_SetSaturation(id, sat)
	# Set cells to print chemistry when print chemistry is turned on
print_chemistry_mask = np.zeros(nxyz,dtype='int')
print_chemistry_mask[:np.int(nxyz/2)]=1
status = RM_SetPrintChemistryMask(id, print_chemistry_mask)
	# Partitioning of uz solids
#status = RM_SetPartitionUZSolids(id, 0)
	# Demonstation of mapping, no symmetry, only one row of cells
grid2chem = -1*np.ones(nxyz)
grid2chem[:np.int(nxyz/2)] = np.arange(np.int(nxyz/2))
grid2chem[np.int(nxyz/2):] = np.arange(np.int(nxyz/2))

grid2chem=grid2chem.astype(np.int32)
status = RM_CreateMapping(id, grid2chem)
if (status < 0):
    status = RM_DecodeError(id, status)
nchem = RM_GetChemistryCellCount(id)




status = RM_SetPrintChemistryOn(id, 0, 1, 0) # workers, initial_phreeqc, utility
	# Set printing of chemistry file
status = RM_LoadDatabase(id, "../database/phreeqc.dat")

	# Demonstration of error handling if ErrorHandlerMode is 0
  # commented out for now (AAE)
	# if (status != IRM_OK)
	# {
	# 	l = RM_GetErrorStringLength(id)
	# 	errstr = (char *) malloc((size_t) (l * sizeof(char) + 1))
	# 	RM_GetErrorString(id, errstr, l+1)
	# 	fprintf(stderr,"Beginning of error string:\n")
	# 	fprintf(stderr,"%s", errstr)
	# 	fprintf(stderr,"End of error string.\n")
	# 	free(errstr)
	# 	errstr = NULL
	# 	RM_Destroy(id)
	# 	exit(1)
	# }
	# Run file to define solutions and reactants for initial conditions, selected output
	# There are three types of IPhreeqc instances in PhreeqcRM
	# Argument 1 refers to the workers for doing reaction calculations for transport
	# Argument 2 refers to the InitialPhreeqc instance for accumulating initial and boundary conditions
	# Argument 3 refers to the Utility instance
status = RM_RunFile(id, 1, 1, 1, "advect.pqi")
	# Clear contents of workers and utility
Str="DELETE -all"
status = RM_RunString(id, 1, 0, 1, Str)	# workers, initial_phreeqc, utility
	# Determine number of components to transport
ncomps = RM_FindComponents(id)

gfw = np.zeros(ncomps,dtype='float64')
status = RM_GetGfw(id, gfw)

#for i in range(ncomps):
# I consider the maximum number of characters for components' name to be 10
components = np.zeros(ncomps,dtype='U10')
for i in range(ncomps):
    status = RM_GetComponent(id, i, components,10)
    print(components[i], "\t", gfw[i], "\n")


ic1 = -1*np.ones(nxyz*7,dtype=np.int32)
ic2 = -1*np.ones(nxyz*7,dtype=np.int32)
f1 = np.ones(nxyz*7,dtype=np.float64)
for i in np.arange(nxyz):
    ic1[i] = 1       # Solution 1
    ic1[i+nxyz] = -1      # Equilibrium phases none
    ic1[i+2*nxyz] = 1       # Exchange 1
    ic1[i+3*nxyz]= -1      # Surface none
    ic1[i+4*nxyz]= -1      # Gas phase none
    ic1[i+5*nxyz]= -1      # Solid solutions none
    ic1[i+6*nxyz]= -1      # Kinetics none
status = RM_InitialPhreeqc2Module(id,ic1, ic2, f1);
t = 0.0 # changed time to time1 to avoid conflicts with julia time function
time_step = 0.0
c = np.zeros(ncomps*nxyz,dtype=np.float64)
status = RM_SetTime(id, t)
status = RM_SetTimeStep(id, time_step)
status = RM_RunCells(id)
status = RM_GetConcentrations(id, c)

nbound = 1
bc1 = np.zeros(nbound,dtype='int32')
bc2 = -1* np.ones(nbound,dtype='int32')
bc_f1 = np.ones(nbound,dtype='float64')
bc_conc = np.zeros(ncomps*nbound,dtype='float64')
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, nbound, bc1, bc2, bc_f1)

#transit loop



	# --------------------------------------------------------------------------
	# Transient loop
	# --------------------------------------------------------------------------
  # this is the beautiful part: I can easily set values after the transport step

nsteps = 10
density = np.ones(nxyz,dtype=np.float64)
volume = np.zeros(nxyz,dtype=np.float64)
pressure = 2*np.ones(nxyz,dtype=np.float64)
temperature = 20*np.ones(nxyz,dtype=np.float64)
sat_calc = np.zeros(nxyz,dtype=np.float64)

status = RM_SetDensity(id, density)
status = RM_SetPressure(id, pressure)
status = RM_SetTemperature(id, temperature)
time_step = 86400
status = RM_SetTimeStep(id, time_step)
for steps in range(nsteps):
    if (steps==nsteps-1):
        print_selected_output_on = 1
        print_chemistry_on = 1
    else:
        print_selected_output_on = 0
        print_chemistry_on = 0
    status = RM_SetSelectedOutputOn(id,print_selected_output_on)
    status = RM_SetPrintChemistryOn(id,print_chemistry_on,0,0)
    status = RM_SetPorosity(id,por)
    status = RM_SetSaturation(id,sat)
    status = RM_SetTemperature(id,temperature)
    status = RM_SetPressure(id,pressure)
    status = RM_SetConcentrations(id,c)
    status = RM_SetTimeStep(id,time_step)
    t = t + time_step
    status = RM_SetTime(id,t)
    strm = 'Begining reaction calculation: '+ str(RM_GetTimeConversion(id))+' days\n'
    strm = strm + '          Time step                         '+ str(RM_GetTime(id) * RM_GetTimeConversion(id))+ ' days \n'
    status = RM_LogMessage(id,strm)
    status = RM_ScreenMessage(id,strm)
    status = RM_RunCells(id)
    status = RM_GetConcentrations(id,c)
    density = np.zeros(nxyz,dtype=np.float64)
    status = RM_GetDensity(id,density)
    volume = np.zeros(nxyz,dtype=np.float64)
    status = RM_GetSolutionVolume(id,volume)
    c1=np.reshape(c.copy(),(nxyz,ncomps),'F')
    if (print_chemistry_on!=0):
        for isel in range(RM_GetSelectedOutputCount(id)):
           #Loop through possible multiple selected output definitions
            n_user = RM_GetNthSelectedOutputUserNumber(id,isel)
            status = RM_SetCurrentSelectedOutputUserNumber(id,n_user)
            print('Selected output sequence number: ',isel, '\n')
            print('Selected output user number:     ',n_user, '\n')
            # Get double array of selected output values

            col = RM_GetSelectedOutputColumnCount(id)
            so = np.zeros(col*nxyz).reshape(nxyz, col)
            status= RM_GetSelectedOutput(id,so)
            for i in range(np.int(RM_GetSelectedOutputRowCount(id)/2)):
                      print('Cell number ', i+1)
                      print('     Density: ',density[i])
                      print('     Volume:  ',volume[i])
                      print('     Components: ')
                      for j in range(ncomps):
                          print('          ',j, ' ', components[j], ': ', c1[i, j])

                      headings = np.zeros(col,dtype='U20')
                      print('     Selected output: ')

                      for j in range(col):
                          status = RM_GetSelectedOutputHeading(id, j,headings[j],20)
                          print('          ',j,' ', headings[j], ': ',so[i,j])


                #status = RM_CloseFiles(id);
                #status = RM_Destroy(id);
