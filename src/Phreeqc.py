from ctypes import *
import numpy as np
import pandas as pd
import sys
libc = CDLL('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/libs/libphreeqcrm-3.5.0.so')
#srcPath = os.path.join(os.path.dirname(__file__),'..')
#sys.path.append('/home/moein/PhreeqcPythonWrapper/PyPhreeqc/src/@IRM_RESULT')
from IRM_RESULT import *
def DECODER(array):
    DecodedArray = np.zeros((1,array.size),dtype='|S20')
    for i in range(array.size):
        DecodedArray[0,i]=array[i].encode()
    return DecodedArray
class PhreeqcRM:
    def __init__(self,nxyz,n_threads):
        self.nxyz = nxyz
        self.n_threads = n_threads
        self.id = libc.RM_Create(nxyz,n_threads)
    def RM_Abort(self,result,err_str):
        return IRM_RESULT(libc.RM_Abort(self.id,result,err_str))

    def RM_CloseFiles(self):
        return libc.RM_CloseFiles(self)

    def RM_Concentrations2Utility(self,c,n,tc,p_atm):
        return libc.RM_Concentrations2Utility(self.id,c,n,tc,p_atm)

    def RM_CreateMapping(self,grid2chem):
        return libc.RM_CreateMapping(self.id,grid2chem.ctypes)

    def RM_DecodeError(self,e):
        return libc.RM_DecodeError(self.id,e)

    def RM_Destroy(self):
        return libc.RM_Destroy(self.id)

    def RM_DumpModule(self,dump_on,append):
        return libc.RM_DumpModule(self.id,dump_on,append)

    def RM_ErrorMessage(self,errstr):
        return libc.RM_ErrorMessage(self.id,errstr)

    def RM_FindComponents(self):
        return libc.RM_FindComponents(self.id)

    def RM_GetBackwardMapping(self,n,list,size):
        return libc.RM_GetBackwardMapping(self.id,n,list,size)

    def RM_GetChemistryCellCount(self):
        return libc.RM_GetChemistryCellCount(self.id)

    def RM_GetComponent(self,num,chem_name,l):
        String = create_string_buffer(l)
        status = libc.RM_GetComponent(self.id,num,String,l)
        chem_name[num] = String.value.decode()
        return status

    def RM_GetConcentrations(self,c):
        return libc.RM_GetConcentrations(self.id,c.ctypes)

    def RM_GetDensity(self,density):
        return libc.RM_GetDensity(self.id,density.ctypes)

    def RM_GetEndCell(self,ec):
        return libc.RM_GetEndCell(self.id,ec)

    def RM_GetEquilibriumPhaseCount(self):
        return libc.RM_GetEquilibriumPhaseCount(self.id)
    def RM_GetEquilibriumPhaseName(self,num,name,l1):
        return libc.RM_GetEquilibriumPhaseName(self.id,num,name,l1)
    def RM_GetErrorString(self,errstr,l):
        return libc.RM_GetErrorString(self.id,errstr,l)
    def RM_GetErrorStringLength(self):
        return libc.RM_GetErrorStringLength(self.id)
    def RM_GetExchangeName(self,num,name,l1):
        return libc.RM_GetExchangeName(self.id,num,name,l1)

    def RM_GetExchangeSpeciesCount(self):
        return libc.RM_GetExchangeSpeciesCount(self.id)

    def RM_GetExchangeSpeciesName(self,num,name,l1):
        return libc.RM_GetExchangeSpeciesName(self.id,num,name,l1)

    def RM_GetFilePrefix(self,prefix,l):
        return libc.RM_GetFilePrefix(self.id,prefix.encode(),l)
    def RM_GetGasComponentsCount(self):
        return libc.RM_GetGasComponentsCount(self.id)

    def RM_GetGasComponentsName(self,nun,name,l1):
        return libc.RM_GetGasComponentsName(self.id,nun,name,l1)

    def RM_GetGfw(self,gfw):
        return libc.RM_GetGfw(self.id,gfw.ctypes)

    def RM_GetGridCellCount(self):
        return libc.RM_GetGridCellCount(self.id)

    def RM_GetIPhreeqcId(self,i):
        return libc.RM_GetIPhreeqcId(self.id,i)

    def RM_GetKineticReactionsCount(self):
        return libc.RM_GetKineticReactionsCount(self.id)

    def RM_GetKineticReactionsName(self,num,name,l1):
        return libc.RM_GetKineticReactionsName(self.id,num,name,l1)

    def RM_GetMpiMyself(self):
        return libc.RM_GetMpiMyself(self.id)

    def RM_GetMpiTasks(self):
        return libc.RM_GetMpiTasks(self.id)

    def RM_GetNthSelectedOutputUserNumber(self,n):
        return libc.RM_GetNthSelectedOutputUserNumber(self.id,n)

    def RM_GetSaturation(self,sat_calc):
        return libc.M_GetSaturation(self.id,sat_calc)

    def RM_GetSelectedOutput(self,so):
        return libc.RM_GetSelectedOutput(self.id,so.ctypes)
    def RM_GetNthSelectedOutputColumnCount(self):
        return libc.RM_GetNthSelectedOutputColumnCount(self.id)

    def RM_GetSelectedOutputCount(self):
        return libc.RM_GetSelectedOutputCount(self.id)

    def RM_GetSelectedOutputHeading(self,num,headings,l):
        String = create_string_buffer(l)
        status = libc.RM_GetSelectedOutputHeading(self.id,icol,String,l)
        headings[num] = String.value.decode()
        return status
    def RM_GetSelectedOutputColumnCount(self):
        return libc.RM_GetSelectedOutputColumnCount(self.id)

    def RM_GetSelectedOutputRowCount(self):
        return libc.RM_GetSelectedOutputRowCount(self.id)

    def RM_GetSICount(self):
        return libc.RM_GetSICount(self.id)

    def RM_GetSIName(self,num,name,l1):
        return libc.RM_GetSIName(self.id,num,name,l1)

    def RM_GetSolidSolutionComponentsCount(self):
        return libc.RM_GetSolidSolutionComponentsCount(self.id)

    def RM_GetSolidSolutionComponentsName(self,num,name,l1):
        return libc.RM_GetSolidSolutionComponentsName(self.id,num,name,l1)

    def RM_GetSolidSolutionName(self,num,name,l1):
        return libc.RM_GetSolidSolutionName(self.id,num,name,l1)

    def RM_GetSolutionVolume(self,vol):
        return libc.RM_GetSolutionVolume(self.id,vol.ctypes)

    def RM_GetSpeciesConcentrations(self,species_conc):
        return libc.RM_GetSpeciesConcentrations(self.id,species_conc)

    def RM_GetSpeciesCount(self):
        return libc.RM_GetSpeciesCount(self.id)

    def RM_GetSpeciesD25(self,diffc):
        return libc.RM_GetSpeciesD25(self.id,diffc)

    def RM_GetSpeciesLog10Gammas(self,species_log10gammas):
        return libc.RM_GetSpeciesLog10Gammas(self.id,species_log10gammas)
    def RM_GetSpeciesName(self,i,name,length):
        return libc.RM_GetSpeciesName(self.id,i,name,length)
    def RM_GetSpeciesSaveOn(self):
        return libc.RM_GetSpeciesSaveOn(self.id)
    def RM_GetSpeciesZ(self,Z):
        return libc.RM_GetSpeciesZ(self.id,Z)
    def RM_GetStartCell(self,sc):
        return libc.RM_GetStartCell(self.id,sc)
    def RM_GetSurfaceName(self,num,name,l1):
        return libc.RM_GetSurfaceName(self.id,num,name,l1)
    def RM_GetSurfaceType(self,num,name,l1):
        return libc.RM_GetSurfaceType(self.id,num,name,l1)
    def RM_GetThreadCount(self):
        return libc.RM_GetThreadCount(self.id)
    def RM_GetTime(self):
        libc.RM_GetTime.restype = c_double
        return libc.RM_GetTime(self.id)
    def RM_GetTimeConversion(self):
        libc.RM_GetTimeConversion.restype = c_double
        return libc.RM_GetTimeConversion(self.id)
    def RM_GetTimeStep(self):
        libc.RM_GetTimeStep.restype = c_double
        return libc.RM_GetTimeStep(self.id)
    def RM_InitialPhreeqc2Module(self,ic1,ic2,f1):
        return libc.RM_InitialPhreeqc2Module(self.id,ic1.ctypes,ic2.ctypes,f1.ctypes)

    def RM_InitialPhreeqc2Concentrations(self,c,n_boundary,boundary_solution1,boundary_solution2,fraction1):
    	return libc.RM_InitialPhreeqc2Concentrations(self.id,c.ctypes,n_boundary,boundary_solution1.ctypes,boundary_solution2.ctypes,fraction1.ctypes)

    def RM_InitialPhreeqc2SpeciesConcentrations(self,species_c,n_boundary,boundary_solution1,boundary_solution2,fraction1):
        return libc.RM_InitialPhreeqc2SpeciesConcentrations(self.id,species_c.ctypes,n_boundary.ctypes,boundary_solution1.ctypes,boundary_solution2.ctypes,fraction1.ctypes)
    def RM_InitialPhreeqcCell2Module(self,n,module_numbers,dim_module_numbers):
        return libc.RM_InitialPhreeqcCell2Module(self.id,n,module_numbers,dim_module_numbers)

    def RM_LoadDatabase(self,db_name):
        return libc.RM_LoadDatabase(self.id,db_name.encode())
    def RM_LogMessage(self,str):
        return libc.RM_LogMessage(self.id,str.encode())
    def RM_MpiWorker(self):
        return libc.RM_MpiWorker(self.id)
    def RM_MpiWorkerBreak(self):
        return libc.RM_MpiWorkerBreak(self.id)
    def RM_OpenFiles(self):
        return libc.RM_OpenFiles(self.id)
    def RM_OutputMessage(self,str):
        return libc.RM_OutputMessage(self.id,str.encode())
    def RM_RunCells(self):
        return libc.RM_RunCells(self.id)
    def RM_RunFile(self,workers,initial_phreeqc,utility,chem_name):
        return libc.RM_RunFile(self.id,workers,initial_phreeqc,utility,chem_name.encode())
    def RM_RunString(self,workers,initial_phreeqc,utility,input_string):
        return libc.RM_RunString(self.id,workers,initial_phreeqc,utility,input_string.encode())
    def RM_ScreenMessage(self,str):
        return libc.RM_ScreenMessage(self.id,str.encode())
    def RM_SetComponentH2O(self,tf):
        return libc.RM_SetComponentH2O(self.id,tf)
    def RM_SetConcentrations(self,c):
        return libc.RM_SetConcentrations(self.id,c.ctypes)
    def RM_SetCurrentSelectedOutputUserNumber(self,n_user):
        return libc.RM_SetCurrentSelectedOutputUserNumber(self.id,n_user)

    def RM_SetDensity(self,density):
        return libc.RM_SetDensity (self.id, density.ctypes)
    def RM_SetDumpFileName(self,dump_name):
        return libc.RM_SetDumpFileName(self.id,dump_name)
    def RM_SetErrorHandlerMode(self,mode):
        return libc.RM_SetErrorHandlerMode(self.id,mode)
    def RM_SetFilePrefix(self,prefix):
        return libc.RM_SetFilePrefix(self.id,prefix.encode())
    def RM_SetMpiWorkerCallbackCookie(self,cookie):
        return libc.RM_SetMpiWorkerCallbackCookie(self.id,cookie)
    def RM_SetPartitionUZSolids(self,tf):
        return libc.RM_SetPartitionUZSolids(self.id,tf)
    def RM_SetPorosity(self, por):
        return libc.RM_SetPorosity(self.id, por.ctypes)
    def RM_SetPressure(self, p):
        return libc.RM_SetPressure(self.id,p.ctypes)
    def RM_SetPrintChemistryMask(self, cell_mask):
        return libc.RM_SetPrintChemistryMask(self.id, cell_mask.ctypes)
    def RM_SetPrintChemistryOn(self,workers,initial_phreeqc,utility):
        return libc.RM_SetPrintChemistryOn(self.id,workers,initial_phreeqc,utility)
    def RM_SetRebalanceByCell(self,method):
        return libc.RM_SetRebalanceByCell(self.id,method)
    def RM_SetRebalanceFraction(self,f):
        return libc.RM_SetRebalanceFraction(self.id,c_double(f))
    def RM_SetRepresentativeVolume(self,rv):
        return libc.RM_SetRepresentativeVolume(self.id,rv.ctypes)
    def RM_SetSaturation(self,sat):
        return libc.RM_SetSaturation(self.id,sat.ctypes)
    def RM_SetScreenOn(self,tf):
        return libc.RM_SetScreenOn(self.id,tf)
    def RM_SetSelectedOutputOn(self,selected_output):
        return libc.RM_SetSelectedOutputOn(self.id,selected_output)
    def RM_SetSpeciesSaveOn(self,save_on):
        return libc.RM_SetSpeciesSaveOn(self.id,save_on)
    def RM_SetTemperature(self,t):
        return libc.RM_SetTemperature(self.id,t.ctypes)
    def RM_SetTime(self,time):
        return libc.RM_SetTime(self.id,c_double(time))
    def RM_SetTimeConversion(self,conv_factor):
        return libc.RM_SetTimeConversion(self.id,c_double(conv_factor))
    def RM_SetTimeStep(self,time_step):
        return libc.RM_SetTimeStep(self.id,c_double(time_step))

    def RM_SetUnitsExchange(self,option):
        return libc.RM_SetUnitsExchange(self.id,option)

    def RM_SetUnitsGasPhase(self,option):
        libc.RM_SetUnitsGasPhase(self.id,option)

    def RM_SetUnitsKinetics(self,option):
        libc.RM_SetUnitsKinetics(self.id,option)

    def RM_SetUnitsPPassemblage(self,option):
        return libc.RM_SetUnitsPPassemblage(self.id,option)

    def RM_SetUnitsSolution(self,option):
        return libc.RM_SetUnitsSolution(self.id,option)

    def RM_SetUnitsSSassemblage(self,option):
        return libc.RM_SetUnitsSSassemblage(self.id,option)
    def RM_SetUnitsSurface(self,option):
        return libc.RM_SetUnitsSurface(self.id,option)

    def RM_SpeciesConcentrations2Module(self,species_conc):
        return libc.RM_SpeciesConcentrations2Module(self.id,species_conc)

    def RM_UseSolutionDensityVolume(self,tf):
        return libc.RM_UseSolutionDensityVolume(self.id,tf)

    def RM_WarningMessage(self,warn_str):
        return libc.RM_WarningMessage(self.id,warn_str)
    def RM_GetComponentCount(self):
        return libc.RM_GetComponentCount(self.id)

    def returnSelectedOutput(self):
        selectedOutput = []
        for isel in range(libc.RM_GetSelectedOutputCount(self.id)):
           #Loop through possible multiple selected output definitions
            n_user = libc.RM_GetNthSelectedOutputUserNumber(self.id,isel)
            status = libc.RM_SetCurrentSelectedOutputUserNumber(self.id,n_user)
            col = libc.RM_GetSelectedOutputColumnCount(self.id)
            so = np.zeros(col*self.nxyz).reshape(self.nxyz, col)
            status= libc.RM_GetSelectedOutput(self.id,so.ctypes)
            headings = np.zeros(col,dtype='U20')

            for j in range(col):
                String = create_string_buffer(20)
                status = libc.RM_GetSelectedOutputHeading(self.id,j,String,20)
                headings[j] = String.value.decode()
            selectedOutput.append(pd.DataFrame(so,columns=headings))
            return selectedOutput
