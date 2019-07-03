import time
T = time.time()
from fipy import Variable, FaceVariable, CellVariable, Grid1D, DiffusionTerm,UpwindConvectionTerm,ConvectionTerm, TransientTerm, DiffusionTerm, Viewer,ExponentialConvectionTerm
from fipy.tools import numerix
from createPhreeqcInstance import *
import pandas as pd
#import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import sys, os
from math import log
import resource
resource.setrlimit(resource.RLIMIT_NOFILE,(65535,65535))
srcPath = os.path.join(os.path.dirname(__file__),'..')
sys.path.append('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/libs')
sys.path.append('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/src/@PhreeqcRM')
sys.path.append('/home/moein/PhreeqcPythonWrapper/Fipy/example/PyPhreeqc/src/@IRM_RESULT')
from PhreeqcRM import *
import matplotlib.pyplot as plt

def extractFromCellVariable(C):
    c = np.zeros(len(C)*C[0].shape[0],dtype='float64')
    nxyz=C[0].shape[0]
    for i in range(len(C)):
        c[i*nxyz:(i+1)*nxyz] = C[i].numericValue
    return c

def putIntoCellVariable(c,C):
    nxyz=C[0].shape[0]
    for i in range(len(C)):
        numerix.put(C[i],np.arange(nxyz),c[i*nxyz:(i+1)*nxyz])
    return C
def yieldCoefficientCalculator(selectedOutput):
    LnQcat = selectedOutput[0]['la_HS-']+2*selectedOutput[0]['la_HCO3-']-selectedOutput[0]['la_SO4-2']-selectedOutput[0]['la_Acetate-']
    LnQana = 0.05*selectedOutput[0]['la_HCO3-']+0.4*selectedOutput[0]['la_H2O']-0.275*selectedOutput[0]['la_H+']-0.2*selectedOutput[0]['la_NH4+']-0.525*selectedOutput[0]['la_Acetate-']
    dgCatSt = -63.9
    dgAnSt = 18.5
    R = 8.314*1e-3
    T = 298.15
    alfa = -0.004
    beta = -0.0694
    v = 0.525
    mm = 0.9306
    dgCat = dgCatSt + R*T*LnQcat
    dgAn = dgAnSt + R*T*LnQana
    #print(dgCat)
    #print(dgAn)
    Y = ((alfa*(dgCatSt**2))-(beta*dgCatSt*dgCat))/((alfa*v*(dgCatSt**2))-(dgCatSt*(beta*v*dgCat+alfa*dgAnSt+dgAn))+(mm*dgCat*dgAnSt))
    dgMet = ((-Y*dgCat*v)+Y*dgAn+dgCat)/Y
    Y[dgCat>-10]=1e9
    if np.min(dgMet) >= 0:
        raise 'metabolism gibbs free energy is not negative anymore!'
    elif np.min(Y) <= 0 :
        raise 'Y is negative'
    return Y
def runReactions(c,time_step,t,phreeqc_rm):
    status = phreeqc_rm.RM_SetConcentrations(c)
    status = phreeqc_rm.RM_SetTimeStep(time_step)
    t = t + time_step
    status = phreeqc_rm.RM_SetTime(t)
    status = phreeqc_rm.RM_RunCells()
    status = phreeqc_rm.RM_GetConcentrations(c)
    status = phreeqc_rm.RM_CloseFiles();
    status = phreeqc_rm.RM_Destroy();
    return c,t
def updatePhreeqcInputFile(Y,filename):
    originalFile = open(filename,'r')
    updatedFile = open("updatedInputFile.phr",'w')
    index = 0
    for line in originalFile:
        if ('KINETICS' in line) and index == 0:
            index=1
            for j in range(100):
                updatedFile.write('KINETICS ' + str(j+1) + '\n')
                updatedFile.write('SO4-2\n')
                updatedFile.write('-formula Acetate -1/'+str(Y[j])+' O 2/'+str(Y[j])+' H 3/'+str(Y[j]) +' C 2/'+str(Y[j])+'\n')
                updatedFile.write('Acetate-\n')
                updatedFile.write('-formula Acetate -0.525 N -0.2 H -0.225 C 0.05 O 0.55 Biomass 1 \n')
                updatedFile.write('-cvode true\n')
        elif index == 1 and '###' in line:
            index = 2
        elif index == 2:
            updatedFile.write(line)
        elif index == 0:
            updatedFile.write(line)
    updatedFile.close()
    originalFile.close()
def checkValidityOfReactiveTransportStep(c,c_old,t,t_old,C,C_old,time_step,tol):
    maxConcentration = np.maximum(c,c_old)
    index = maxConcentration>1e-9
    relError=np.abs((c[index]-c_old[index])/maxConcentration[index])
    if np.max(relError)>tol:
        C=list(C_old)
        t=t_old
        time_step=time_step/2
        print('the step was not accepted, comming back')
        index=0
    else:
        print('the step was accepted lets move to the next step')
        time_step=time_step*1.1
        C = putIntoCellVariable(c,C)
        index=1

    return C, t, time_step,index
def plot(filename1,filename2):
    withGEDYM=pd.read_csv(filename1+'.csv')
    withoutGEDYM=pd.read_csv(filename2+'.csv')
    withGEDYM.rename(columns={'Unnamed: 0':'index','0':'time(day)','1':'Acetate Concentration(M)'},inplace=True)
    withoutGEDYM.rename(columns={'Unnamed: 0':'index','0':'time(day)','1':'Acetate Concentration(M)'},inplace=True)
    withGEDYM.set_index('index',inplace=True)
    withoutGEDYM.set_index('index',inplace=True)
    plt.figure(1)
    ax1=plt.scatter(x='time(day)',y='Acetate Concentration(M)',data=withGEDYM,c='b')
    plt.xlim(withGEDYM['time(day)'].min(),withGEDYM['time(day)'].max())
    plt.ylim(withGEDYM['Acetate Concentration(M)'].min()-0.001,withGEDYM['Acetate Concentration(M)'].max()+0.001)
    ax2=plt.scatter(x='time(day)',y='Acetate Concentration(M)',data=withoutGEDYM,c='r')
    plt.xlabel('time(day)')
    plt.ylabel('Acetate Concentration(M)')
    plt.figure(2)
    a=np.abs(withGEDYM['Acetate Concentration(M)']-withoutGEDYM['Acetate Concentration(M)'])*100/withGEDYM['Acetate Concentration(M)']
    b=ithGEDYM['Acetate Concentration(M)']
    return a,b
#Acetate-
 # -formula Acetate -0.525 N -0.2 H -0.225 C 0.05 O 0.55 Biomass 1
  #-m 0







#transit loop



	# --------------------------------------------------------------------------
	# Transient loop
	# --------------------------------------------------------------------------
  # this is the beautiful part: I can easily set values after the transport step

maxTimeStep=120
clConcentration =  np.zeros(maxTimeStep,dtype='float64')
zaman =  np.zeros(maxTimeStep,dtype='float64')

bc_conc,mesh,ncomps,components,c,nxyz,phreeqc_rm_toCalculateActivityCoef = createPhreeqcInstance("onlySolutionData.phr",-1)

nSpecies = phreeqc_rm_toCalculateActivityCoef.RM_GetSpeciesCount()

#==================================================
#lets set up the transport model
C=[]
for i in range(ncomps):
    C.append(CellVariable(name=components[i].encode(),mesh=mesh,value=c[i*nxyz:(i+1)*nxyz]))
    C[i].constrain(bc_conc[i],mesh.facesLeft)
    C[i].faceGrad.constrain(0,mesh.facesRight)
if __name__ == '__main__':
    viewer = Viewer(vars=(C[4]),datamin=0.,datamax=1e-4)
    viewer.plot()
#===================================================
D = 0.0 # assume there is no diffusion except numerical diffusion hahahaha
u = (1.5/432000,) # m/s equivalent to 1 dt/day
eqX = []
for i in range(ncomps):
    if (i!=4):
        eqX.append(TransientTerm() == DiffusionTerm(coeff=D)-UpwindConvectionTerm(coeff=u))
    else:
        # I am assuming that microbes are attached
        #eqX.append(TransientTerm() == DiffusionTerm(coeff=D))
        eqX.append(TransientTerm() == DiffusionTerm(coeff=D)-UpwindConvectionTerm(coeff=u))
#timeStepDuration = 0.9 * dx**2 / (2 * D) # the maximum time step
time_step = 432000.0
t=0
C_old = list(C)
steps=-1
t_final = 120 * 432000.0
tol=0.1
trial = 0
H2SProductionData=np.zeros((maxTimeStep,nxyz),dtype='float')
while (t<t_final):
    trial +=1
    #=======================================================================
    # Transport Seteps
    #cBeforTransport=c.copy()
    for i in range(ncomps):
        eqX[i].solve(var=C[i],dt=time_step)
    c = extractFromCellVariable(C)
    #cAfterTransport=c.copy()
    #=========================================================================
    # Reaction steps
    # 1) calculate yield yieldCoefficient
    # export concentrations to the yieldCoefficientCaluclator instances
    status = phreeqc_rm_toCalculateActivityCoef.RM_SetConcentrations(c)
    c_old = c.copy()
    status = phreeqc_rm_toCalculateActivityCoef.RM_RunCells()
    # Recieve activities and calculate yield coefficients
    selectedOutput=phreeqc_rm_toCalculateActivityCoef.RM_ReturnSelectedOutput
    withGEDYM = False
    if withGEDYM == True :

        Y = yieldCoefficientCalculator(selectedOutput)
    # Based on activities create the suitable input file for the kinetic calculations
        updatePhreeqcInputFile(Y,"with NRB.phr")
    # now create the instance for the kinetic reaction calculations
        bc_conc_fake,mesh_fake,ncomps_fake,components_fake,c_fake,nxyz_fake,phreeqc_rm_KineticInstance = createPhreeqcInstance("updatedInputFile.phr",1)
    else:
        bc_conc_fake,mesh_fake,ncomps_fake,components_fake,c_fake,nxyz_fake,phreeqc_rm_KineticInstance = createPhreeqcInstance("with NRB.phr",1)
    # now its time to set the concentrations and time and time step and run the reactions for this instancee
    t_old=t
    c, t = runReactions(c, time_step, t, phreeqc_rm_KineticInstance)
    C = putIntoCellVariable(c,C)
    # here we check, if the time step is so long, we make it smaller to make Sure
    # that in each time step reactions are not changing the concentration more than 5%

    print('======================================================')
    print(t)
    print('simulation is ' + np.str(t/t_final) +' completed')

    steps += 1
    zaman[steps]=t/(24*3600)
    H2SProductionData[steps]=C[3].numericValue
    plot('withGEDYM','withoutGEDYM')

    #viewer.plot()
    # in this block you have to write a piece of code that measures the activity coefficient
    #selectedOutput = phreeqc_rm_toCalculateActivityCoef.returnSelectedOutput()


    #c,t,phreeqc_rm_toCalculateActivityCoef=runReactions(c,time_step,t,phreeqc_rm_toCalculateActivityCoef)
    #C = putIntoCellVariable(c,C)
    #
    #if __name__ == '__main__':
    #    viewer.plot()
    #clConcentration[steps]=C[4].numericValue[99]
