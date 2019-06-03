from ctypes import *
import numpy as np
libc = CDLL('../libs/libphreeqcrm.so')

def RM_Abort(id,result,err_str):
    return (libc.RM_Abort(id,result,err_str))

def RM_CloseFiles(id):
    return libc.RM_CloseFiles(id)

def RM_Concentrations2Utility(id,c,n,tc,p_atm):
    return libc.RM_Concentrations2Utility(id,c,n,tc,p_atm)


def RM_Create(nxyz,n_threads):
#    if not(is_int(nxyz) or is_int(n_threads)):
#        raise Exception('nxyz and n_threads should be integers')
    return libc.RM_Create(nxyz,n_threads)

def RM_CreateMapping(id,grid2chem):
    return libc.RM_CreateMapping(id,grid2chem.ctypes)

def RM_DecodeError(id,e):
    return libc.RM_DecodeError(id,e)

def RM_Destroy(id):
    return libc.RM_Destroy(id)

def RM_DumpModule(id,dump_on,append):
    return libc.RM_DumpModule(id,dump_on,append)

def RM_ErrorMessage(id,errstr):
    return libc.RM_ErrorMessage(id,errstr)

def RM_FindComponents(id):
    return libc.RM_FindComponents(id)

def RM_GetBackwardMapping(id,n,list,size):
    return libc.RM_GetBackwardMapping(id,n,list,size)

def RM_GetChemistryCellCount(id):
    return libc.RM_GetChemistryCellCount(id)

def RM_GetComponent(id,num,chem_name,l):
    String = create_string_buffer(l)
    status = libc.RM_GetComponent(id,num,String,l)
    chem_name[num] = String.value.decode()
    return status

def RM_GetConcentrations(id,c):
    return libc.RM_GetConcentrations(id,c.ctypes)

def RM_GetDensity(id,density):
    return libc.RM_GetDensity(id,density.ctypes)

def RM_GetEndCell(id,ec):
    return libc.RM_GetEndCell(id,ec)

def RM_GetEquilibriumPhaseCount(id):
    return libc.RM_GetEquilibriumPhaseCount(id)
def RM_GetEquilibriumPhaseName(id,num,name,l1):
    return libc.RM_GetEquilibriumPhaseName(id,num,name,l1)
def RM_GetErrorString(id,errstr,l):
    return libc.RM_GetErrorString(id,errstr,l)
def RM_GetErrorStringLength(id):
    return libc.RM_GetErrorStringLength(id)
def RM_GetExchangeName(id,num,name,l1):
    return libc.RM_GetExchangeName(id,num,name,l1)

def RM_GetExchangeSpeciesCount(id):
    return libc.RM_GetExchangeSpeciesCount(id)

def RM_GetExchangeSpeciesName(id,num,name,l1):
    return libc.RM_GetExchangeSpeciesName(id,num,name,l1)

def RM_GetFilePrefix(id,prefix,l):
    return libc.RM_GetFilePrefix(id,prefix.encode(),l)
def RM_GetGasComponentsCount(id):
    return libc.RM_GetGasComponentsCount(id)

def RM_GetGasComponentsName(id,nun,name,l1):
    return libc.RM_GetGasComponentsName(id,nun,name,l1)

def RM_GetGfw(nxyz,gfw):
    return libc.RM_GetGfw(nxyz,gfw.ctypes)

def RM_GetGridCellCount(id):
    return libc.RM_GetGridCellCount(id)

def RM_GetIPhreeqcId(id,i):
    return libc.RM_GetIPhreeqcId(id,i)

def RM_GetKineticReactionsCount(id):
    return libc.RM_GetKineticReactionsCount(id)

def RM_GetKineticReactionsName(id,num,name,l1):
    return libc.RM_GetKineticReactionsName(id,num,name,l1)

def RM_GetMpiMyself(id):
    return libc.RM_GetMpiMyself(id)

def RM_GetMpiTasks(id):
    return libc.RM_GetMpiTasks(id)

def RM_GetNthSelectedOutputUserNumber(id,n):
    return libc.RM_GetNthSelectedOutputUserNumber(id,n)

def RM_GetSaturation(id,sat_calc):
    return libc.M_GetSaturation(id,sat_calc)

def RM_GetSelectedOutput(id,so):
    return libc.RM_GetSelectedOutput(id,so.ctypes)
def RM_GetNthSelectedOutputColumnCount(id):
    return libc.RM_GetNthSelectedOutputColumnCount(id)

def RM_GetSelectedOutputCount(id):
    return libc.RM_GetSelectedOutputCount(id)

def RM_GetSelectedOutputHeading(id,icol,heading,length):
    return libc.RM_GetSelectedOutputHeading(id,icol,heading,length)
def RM_GetSelectedOutputColumnCount(id):
    return libc.RM_GetSelectedOutputColumnCount(id)

def RM_GetSelectedOutputRowCount(id):
    return libc.RM_GetSelectedOutputRowCount(id)

def RM_GetSICount(id):
    return libc.RM_GetSICount(id)

def RM_GetSIName(id,num,name,l1):
    return libc.RM_GetSIName(id,num,name,l1)

def RM_GetSolidSolutionComponentsCount(id):
    return libc.RM_GetSolidSolutionComponentsCount(id)

def RM_GetSolidSolutionComponentsName(id,num,name,l1):
    return libc.RM_GetSolidSolutionComponentsName(id,num,name,l1)

def RM_GetSolidSolutionName(id,num,name,l1):
    return libc.RM_GetSolidSolutionName(id,num,name,l1)

def RM_GetSolutionVolume(id,vol):
    return libc.RM_GetSolutionVolume(id,vol.ctypes)

def RM_GetSpeciesConcentrations(id,species_conc):
    return libc.RM_GetSpeciesConcentrations(id,species_conc)

def RM_GetSpeciesCount(id):
    return libc.RM_GetSpeciesCount(id)

def RM_GetSpeciesD25(id,diffc):
    return libc.RM_GetSpeciesD25(id,diffc)

def RM_GetSpeciesLog10Gammas(id,species_log10gammas):
    return libc.RM_GetSpeciesLog10Gammas(id,species_log10gammas)
def RM_GetSpeciesName(id,i,name,length):
    return libc.RM_GetSpeciesName(id,i,name,length)
def RM_GetSpeciesSaveOn(id):
    return libc.RM_GetSpeciesSaveOn(id)
def RM_GetSpeciesZ(id,Z):
    return libc.RM_GetSpeciesZ(id,Z)
def RM_GetStartCell(id,sc):
    return libc.RM_GetStartCell(id,sc)
def RM_GetSurfaceName(id,num,name,l1):
    return libc.RM_GetSurfaceName(id,num,name,l1)
def RM_GetSurfaceType(id,num,name,l1):
    return libc.RM_GetSurfaceType(id,num,name,l1)
def RM_GetThreadCount(id):
    return libc.RM_GetThreadCount(id)
def RM_GetTime(id):
    libc.RM_GetTime.restype = c_double
    return libc.RM_GetTime(id)
def RM_GetTimeConversion(id):
    libc.RM_GetTimeConversion.restype = c_double
    return libc.RM_GetTimeConversion(id)
def RM_GetTimeStep(id):
    libc.RM_GetTimeStep.restype = c_double
    return libc.RM_GetTimeStep(id)
def RM_InitialPhreeqc2Module(id,ic1,ic2,f1):
    return libc.RM_InitialPhreeqc2Module(id,ic1.ctypes,ic2.ctypes,f1.ctypes)

def RM_InitialPhreeqc2Concentrations(id,c,n_boundary,boundary_solution1,boundary_solution2,fraction1):
	return libc.RM_InitialPhreeqc2Concentrations(id,c.ctypes,n_boundary,boundary_solution1.ctypes,boundary_solution2.ctypes,fraction1.ctypes)

def RM_InitialPhreeqc2SpeciesConcentrations(id,species_c,n_boundary,boundary_solution1,boundary_solution2,fraction1):
    return libc.RM_InitialPhreeqc2SpeciesConcentrations(id,species_c.ctypes,n_boundary.ctypes,boundary_solution1.ctypes,boundary_solution2.ctypes,fraction1.ctypes)
def RM_InitialPhreeqcCell2Module(id,n,module_numbers,dim_module_numbers):
    return libc.RM_InitialPhreeqcCell2Module(id,n,module_numbers,dim_module_numbers)

def RM_LoadDatabase(id,db_name):
    return libc.RM_LoadDatabase(id,db_name.encode())
def RM_LogMessage(id,str):
    return libc.RM_LogMessage(id,str.encode())
def RM_MpiWorker(id):
    return libc.RM_MpiWorker(id)
def RM_MpiWorkerBreak(id):
    return libc.RM_MpiWorkerBreak(id)
def RM_OpenFiles(id):
    return libc.RM_OpenFiles(id)
def RM_OutputMessage(id,str):
    return libc.RM_OutputMessage(id,str.encode())
def RM_RunCells(id):
    return libc.RM_RunCells(id)
def RM_RunFile(id,workers,initial_phreeqc,utility,chem_name):
    return libc.RM_RunFile(id,workers,initial_phreeqc,utility,chem_name.encode())
def RM_RunString(id,workers,initial_phreeqc,utility,input_string):
    return libc.RM_RunString(id,workers,initial_phreeqc,utility,input_string.encode())
def RM_ScreenMessage(id,str):
    return libc.RM_ScreenMessage(id,str.encode())
def RM_SetComponentH2O(id,tf):
    return libc.RM_SetComponentH2O(id,tf)
def RM_SetConcentrations(id,c):
    return libc.RM_SetConcentrations(id,c.ctypes)
def RM_SetCurrentSelectedOutputUserNumber(id,n_user):
    return libc.RM_SetCurrentSelectedOutputUserNumber(id,n_user)

def RM_SetDensity(id,density):
    return libc.RM_SetDensity (id, density.ctypes)
def RM_SetDumpFileName(id,dump_name):
    return libc.RM_SetDumpFileName(id,dump_name)
def RM_SetErrorHandlerMode(id,mode):
    return libc.RM_SetErrorHandlerMode(id,mode)
def RM_SetFilePrefix(id,prefix):
    return libc.RM_SetFilePrefix(id,prefix.encode())
def RM_SetMpiWorkerCallbackCookie(id,cookie):
    return libc.RM_SetMpiWorkerCallbackCookie(id,cookie)
def RM_SetPartitionUZSolids(id,tf):
    return libc.RM_SetPartitionUZSolids(id,tf)
def RM_SetPorosity(id, por):
    return libc.RM_SetPorosity(id, por.ctypes)
def RM_SetPressure(id, p):
    return libc.RM_SetPressure(id,p.ctypes)
def RM_SetPrintChemistryMask(id, cell_mask):
    return libc.RM_SetPrintChemistryMask(id, cell_mask.ctypes)
def RM_SetPrintChemistryOn(id,workers,initial_phreeqc,utility):
    return libc.RM_SetPrintChemistryOn(id,workers,initial_phreeqc,utility)
def RM_SetRebalanceByCell(id,method):
    return libc.RM_SetRebalanceByCell(id,method)
def RM_SetRebalanceFraction(id,f):
    return libc.RM_SetRebalanceFraction(id,c_double(f))
def RM_SetRepresentativeVolume(id,rv):
    return libc.RM_SetRepresentativeVolume(id,rv.ctypes)
def RM_SetSaturation(id,sat):
    return libc.RM_SetSaturation(id,sat.ctypes)
def RM_SetScreenOn(id,tf):
    return libc.RM_SetScreenOn(id,tf)
def RM_SetSelectedOutputOn(id,selected_output):
    return libc.RM_SetSelectedOutputOn(id,selected_output)
def RM_SetSpeciesSaveOn(id,save_on):
    return libc.RM_SetSpeciesSaveOn(id,save_on)
def RM_SetTemperature(id,t):
    return libc.RM_SetTemperature(id,t.ctypes)
def RM_SetTime(id,time):
    return libc.RM_SetTime(id,c_double(time))
def RM_SetTimeConversion(id,conv_factor):
    return libc.RM_SetTimeConversion(id,c_double(conv_factor))
def RM_SetTimeStep(id,time_step):
    return libc.RM_SetTimeStep(id,c_double(
    ))

def RM_SetUnitsExchange(id,option):
    return libc.RM_SetUnitsExchange(id,option)

def RM_SetUnitsGasPhase(id,option):
    libc.RM_SetUnitsGasPhase(id,option)

def RM_SetUnitsKinetics(id,option):
    libc.RM_SetUnitsKinetics(id,option)

def RM_SetUnitsPPassemblage(id,option):
    return libc.RM_SetUnitsPPassemblage(id,option)

def RM_SetUnitsSolution(id,option):
    return libc.RM_SetUnitsSolution(id,option)

def RM_SetUnitsSSassemblage(id,option):
    return libc.RM_SetUnitsSSassemblage(id,option)

def RM_SetUnitsSurface(id,option):
    return libc.RM_SetUnitsSurface(id,option)

def RM_SpeciesConcentrations2Module(id,species_conc):
    return libc.RM_SpeciesConcentrations2Module(id,species_conc)

def RM_UseSolutionDensityVolume(id,tf):
    return libc.RM_UseSolutionDensityVolume(id,tf)

def RM_WarningMessage(id,warn_str):
    return libc.RM_WarningMessage(id,warn_str)
def RM_GetComponentCount(id):
    return libc.RM_GetComponentCount(id)
