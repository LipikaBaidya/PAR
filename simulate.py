#!/usr/bin/env python
# coding: utf-8

# In[1]:


from simtk import openmm as mm
from simtk import unit as u
from simtk.openmm import app
import numpy as np
import math as m
import os, sys
import datetime
import json
from sys import platform, stdout
import time


time_st = datetime.datetime.now()
print("Simulation started at \t",)
print(time_st)


#inputfile = json.load(open('input.json'))
inputfile = json.load(open(sys.argv[1]))


unit_len = 0.1*u.nanometer #1 Angstrom


# In[2]:


##### ``````````````` Reading Input parameters```````````````````#####

nom       = inputfile['nom']
nPos      = inputfile['nPos']
nNeg      = inputfile['nNeg']
nMg       = inputfile['nMg']
mPos      = inputfile['mPos']
mNeg      = inputfile['mNeg']
mMg       = inputfile['mMg']
totstep   = inputfile['totstep']
timeint   = inputfile['timeint']
kB        = inputfile['kB'] 
sigma_bb  = inputfile['sigma_bb']
epsilon_l = inputfile['epsilon_l']
T         = inputfile['T'] 
kb_PP     = inputfile['kb_PP']
kb_AS     = inputfile['kb_AS']
kb_PS     = inputfile['kb_PS']
kb_SP     = inputfile['kb_SP']
kb_SS     = inputfile['kb_SS']
r0_PP     = inputfile['r0_PP']
r0_AS     = inputfile['r0_AS']
r0_PS     = inputfile['r0_PS']
r0_SP     = inputfile['r0_SP']
r0_SS     = inputfile['r0_SS']
ka_PPS    = inputfile['ka_PPS']
ka_BSP    = inputfile['ka_BSP']
ka_PSS    = inputfile['ka_PSS']
ka_BSS    = inputfile['ka_BSS']
ka_SSP1   = inputfile['ka_SSP1']
ka_SP1P1  = inputfile['ka_SP1P1']
theta0_PPS   = inputfile['theta0_PPS']
theta0_BSP   = inputfile['theta0_BSP']
theta0_PSS   = inputfile['theta0_PSS']
theta0_BSS   = inputfile['theta0_BSS']
theta0_SSP1  = inputfile['theta0_SSP1']
theta0_SP1P1 = inputfile['theta0_SP1P1']
h_ST         = inputfile['h_ST']
Tm           = inputfile['Tm']
S_ST         = inputfile['S_ST']
krST         = inputfile['krST']
kphi1ST      = inputfile['kphi1ST']
kphi2ST      = inputfile['kphi2ST']
R0ST         = inputfile['R0ST']
phi10ST      = inputfile['phi10ST']
phi20ST      = inputfile['phi20ST']
R0        = inputfile['R0']
salt_conc = inputfile['salt_conc']
epsilon_r = inputfile['epsilon_r']
omega     = inputfile['omega']
cutoffNL  = inputfile['rcutNL']
rcutElec  = inputfile['rcutElec']
epBB      = inputfile['epBB']
epBB_baul = inputfile['epBB_baul']
epBS_baul = inputfile['epBS_baul']
epSS_baul = inputfile['epSS_baul']
box       = inputfile['box']
sigmaPos  = inputfile['sigmaPos']
sigmaNeg  = inputfile['sigmaNeg']
sigmaMg   = inputfile['sigmaMg']
epsPos    = inputfile['epsilonPos']
epsNeg    = inputfile['epsilonNeg']
epsMg     = inputfile['epsilonMg']
epsPI     = inputfile['epsPI']
qPos      = inputfile['qPos']
qNeg      = inputfile['qNeg']
qMg       = inputfile['qMg']

ntot      = nom + nPos + nNeg + nMg


#*************   Assigning radii value for each bead *******************#
sigma = []
mass  = []
for i in range(0,nom):
	if (((i%5) == 0) or ((i%5) == 1)): #P
		sigma_par = 1.89
		mass_par  = 78.9596 * u.amu
	elif ((i%5) == 2) or ((i%5) == 4): #S
		sigma_par = 2.61
		mass_par  = 117.0557 * u.amu
	else:
		sigma_par = 2.52 #B
		mass_par  = 134.0472 * u.amu
	sigma.append(sigma_par)
	mass.append(mass_par)
#for i in range(0,nom):
#	print(i,sigma[i],mass[i])

for i in range(nom,nom+nNeg):
	sigma.append(sigmaNeg)
	mass.append(mNeg *u.amu)
for i in range(nom+nNeg,nom+nNeg+nPos):
	sigma.append(sigmaPos)
	mass.append(mPos *u.amu)
for i in range(nom+nNeg+nPos,ntot):
        sigma.append(sigmaMg)
        mass.append(mMg *u.amu)
#### ```````````` Converting in openmm units ```````````````````####


T         = T * u.kelvin
unit_kb   = mm.KJPerKcal * u.kilojoules_per_mole/unit_len**2
epsilon_l = epsilon_l * mm.KJPerKcal * u.kilojoule_per_mole 
cutoffNL  = cutoffNL * unit_len
box       = box * unit_len 

tau_L     = 1.0 * u.picoseconds 
fric_coef = 1.0/u.picoseconds
timeStep  = 2.0 * u.femtoseconds

### ````````````` Printing Input values ```````````````````````````####

print("Number of particles: protein, positive ion, negative ion and total particles\t")
#print(nom,nPos,nNeg,ntot)



#### `````````````` Information of electrostatic parameters (Debye Huckel) #####

lb_kT    = (332.063709/epsilon_r) * unit_len * mm.KJPerKcal * u.kilojoules_per_mole
num_den  = 2.0 * salt_conc * 6.022 * 0.0001/(unit_len)**3
kappa    = (4.0 * m.pi * lb_kT * num_den/(kB * T))**0.5


# In[3]:


################## Input coordinates for simulations and creating system object #####################################

pdb = app.PDBFile("pos.pdb")
positions = pdb.positions
#print((positions))
system = mm.System()
system.setDefaultPeriodicBoxVectors([box,0,0], [0,box,0], [0,0,box])


for i in range(ntot):
	system.addParticle(mass[i])
# In[4]:


##`````````````````````   Adding Forces  ''''''''''''''''''''''''''''''##

#^^^^^^^^^^^^^^^^^   Bonded FORCE ^^^^^^^^^^^^^^^^^^^^^^^^^^ #

feneBond = mm.HarmonicBondForce()

for i in range(0,int(nom/5)):
	iP1 = 5*i    # P1
	iP2 = 5*i+1  # P2
	iS1 = 5*i+2  # S
	iB  = 5*i+3  # B1
	iS2 = 5*i+4  # S2
	feneBond.addBond(iP1,iP2, r0_PP * unit_len, 2.0 * kb_PP * unit_kb)  # P1-P2
	feneBond.addBond(iP2,iS1, r0_PS * unit_len, 2.0 * kb_PS * unit_kb)  # P1-S1
	feneBond.addBond(iS1,iB, r0_AS * unit_len, 2.0 * kb_AS * unit_kb) # B-S1
	feneBond.addBond(iS1,iS2, r0_SS * unit_len, 2.0 * kb_SS * unit_kb) # S1-S2
	if(i < int(nom/5)-1 ):
		feneBond.addBond(iS2,iS2+1, r0_SP * unit_len, 2.0 * kb_SP * unit_kb) # S2-(P1+1)
system.addForce(feneBond)

#^^^^^^^^^^^^^^^^^   Angle FORCE ^^^^^^^^^^^^^^^^^^^^^^^^^^ #
harmonicAngle = mm.HarmonicAngleForce()
for i in range(0,int(nom/5)):
	iP1 = 5*i    # P1
	iP2 = 5*i+1  # P2
	iS1 = 5*i+2  # S1
	iB  = 5*i+3  # B
	iS2 = 5*i+4  # S2
	iP1_1 = iS2+1 # P1+1
	iP2_1 = iS2+2 # P2+1
	harmonicAngle.addAngle(iP1, iP2, iS1, theta0_PPS*u.radians, 2.0 * ka_PPS*u.kilocalories_per_mole/u.radians**2) # P1-P2-S1
	harmonicAngle.addAngle(iB, iS1, iP2, theta0_BSP*u.radians, 2.0 * ka_BSP*u.kilocalories_per_mole/u.radians**2)  # B-S1-P2
	harmonicAngle.addAngle(iP2, iS1, iS2, theta0_PSS*u.radians, 2.0 * ka_PSS*u.kilocalories_per_mole/u.radians**2)  # P-S1-S2
	harmonicAngle.addAngle(iB, iS1, iS2, theta0_BSS*u.radians, 2.0 * ka_BSS*u.kilocalories_per_mole/u.radians**2)  # B-S1-S2
	if(i < int(nom/5)-1 ):
		harmonicAngle.addAngle(iS1, iS2, iP1_1, theta0_SSP1*u.radians, 2.0 * ka_SSP1*u.kilocalories_per_mole/u.radians**2) # S1-S2-(P1+1)
		harmonicAngle.addAngle(iS2, iP1_1, iP2_1, theta0_SP1P1*u.radians, 2.0 * ka_SP1P1*u.kilocalories_per_mole/u.radians**2)  # S2-(P1+1)-(P2+1)
system.addForce(harmonicAngle)


#^^^^^^^^^^^^^^^ ELECTROSTATIC FORCE (COLUMB) AMONG ALL BEADS (protein-protein, protein-ion, ion-ion) ^^^^^^^^^^^^^^^^^^^^^#
#~~~~~~~~~~~~~~~~  Assigning charge for all beads ~~~~~~~~~~~~~#
q_bead = []
for i in range(nom):
	if ((i%5)==0 or (i%5)==1):
		q_bead.append(-1.0)
	else:
		q_bead.append(0.0)
	#print(i,q_bead[i])

for i in range(nom,nom+nNeg):
    q_bead.append(qNeg)
for i in range(nom+nNeg,nom+nNeg+nPos):
    q_bead.append(qPos)
for i in range(nom+nNeg+nPos,ntot):
    q_bead.append(qMg)
#~~~~~~~~~~~~~~~~~~ Debye-Huckel Interaction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


KELVIN_TO_KT 	= u.AVOGADRO_CONSTANT_NA * u.BOLTZMANN_CONSTANT_kB / u.kilocalorie_per_mole
b_DH 		= 4.38178046 * u.angstrom / u.elementary_charge
T_unitless 	= T * KELVIN_TO_KT
epsilonDH 	= 296.0736276 - 619.2813716 * T_unitless + 531.2826741 * T_unitless**2 - 180.0369914 * T_unitless**3;
l_Bjerrum 	= (332.0637*u.angstroms)/epsilonDH
salt_conc_mM    = salt_conc * 1000.0
kappaDH 	= u.sqrt(4* np.pi * l_Bjerrum * 2 * salt_conc_mM * 6.022e-7 / (T_unitless * u.angstrom**3))
Q_DH 		= b_DH * T_unitless * u.elementary_charge**2 / l_Bjerrum
print("KELVIN_TO_KT:",KELVIN_TO_KT)
print("b_DH:",b_DH)
print("T_unitless:",T_unitless)
print("epsilonDH:",epsilonDH)
print("l_Bjerrum:",l_Bjerrum)
print("salt_conc_mM:",salt_conc_mM)
print("kappaDH:",kappaDH)
print("Q_DH:",Q_DH)
print("scaleDH:",l_Bjerrum * u.kilocalorie_per_mole / u.elementary_charge**2)

b_unitless  = b_DH/(u.angstrom / u.elementary_charge)
b3 = b_unitless**3
Q_unitless = Q_DH/u.elementary_charge
V2 = 4 * np.pi * np.e * b3 * (1. + qMg) * (1 / Q_unitless - 1./ qMg)
V1 = 4 * np.pi * np.e * b3 * (1 + 1) * (1/Q_unitless - 1)
C1V1 = salt_conc_mM * V1 * 6.022e-07
#print("C1V1 = ", C1V1)
Mg_conc_in_mM = nMg /( 6.022e-7 * (box/unit_len)**3)
print("Mg_conc_in_mM:",Mg_conc_in_mM)
#N_Mg = int(args.divalent_concentration * 6.022e-7 * pow(args.box_size, 3))
C2V2 = Mg_conc_in_mM * V2 * 6.022e-07
#print("C2V2 = ", C2V2)

theta1 = 0
if(Mg_conc_in_mM > 0):
    if (qMg == 2):
        theta1 = (np.sqrt(pow(C1V1, 4) + 8 * np.e * C1V1 * C1V1 * C2V2 * (1.0 - Q_unitless)) - C1V1 * C1V1) / (4 * np.e * C2V2)
        Q_DH   = (1.0 - theta1) * u.elementary_charge
print("Phosphate charge after updates = ", Q_DH)

typeDH = []
for i in range(nom):
    typeDH.append(1)
for i in range(nom,ntot):
    typeDH.append(-1)

DHforce = mm.CustomNonbondedForce("0.5 * abs(typeDH1 + typeDH2) * scaleDH * q1 * q2 * exp(-kappaDH*r)/r")
DHforce.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
DHforce.addPerParticleParameter('typeDH')
DHforce.addGlobalParameter("scaleDH", l_Bjerrum * u.kilocalorie_per_mole / u.elementary_charge**2)
DHforce.addGlobalParameter("kappaDH", kappaDH)
DHforce.addPerParticleParameter('q')
DHforce.setCutoffDistance(0.5*box)

for i in range(nom):     
        DHforce.addParticle([typeDH[i],q_bead[i] * Q_DH])

for i in range(nom,ntot):  
        DHforce.addParticle([typeDH[i],q_bead[i] * u.elementary_charge])

for i in range(0,int(nom/5)):
	iP1 = 5*i    # P1
	iP2 = 5*i+1  # P2
	iS1 = 5*i+2  # S1
	iB  = 5*i+3  # B
	iS2 = 5*i+4  # S2
	DHforce.addExclusion(iP1,iP2)
	DHforce.addExclusion(iP2,iS1)
	DHforce.addExclusion(iS1,iB)
	DHforce.addExclusion(iS1,iS2)
	if(i < int(nom/5)-1 ):
		DHforce.addExclusion(iS2,iS2+1)

system.addForce(DHforce)
#~~~~~~~~~~~~~~~~ Mg - P RISM Potential ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
switchMP = [
  # P               B/S             Mg 
    0,              0,              1,      # P  
    0,              0,              0,      # B/S
    1,              0,              0]      # Mg2+

TypeMP = []
for i in range(0,nom):
    if ((i%5) == 0 or (i%5) == 1): #P
            TypeMP.append(0)
    else:
            TypeMP.append(1) #B/S
for i in range(nom,ntot):
    TypeMP.append(2)

MP_potential = mm.CustomNonbondedForce("switchMP(typeMP1,typeMP2) * (W(r) + ((D_H - W(r)) * exp(-a^2/r^2))); D_H=scaleDH * q1 * q2 * exp(-kappaDH*r)/r")
MP_potential.addTabulatedFunction('switchMP', mm.Discrete2DFunction(3, 3, switchMP))
MP_potential.addPerParticleParameter('typeMP')
MP_potential.addGlobalParameter("scaleDH", l_Bjerrum * u.kilocalorie_per_mole / u.elementary_charge**2)
MP_potential.addGlobalParameter("kappaDH", kappaDH)
MP_potential.addGlobalParameter("a", 5.0 * u.angstroms)
MP_potential.addPerParticleParameter('q')
MP_potential.setCutoffDistance(0.5*box)
MP_potential.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)

with open('Mg_P_potential.dat', 'r') as file:
    tabulated_data = [float(line.strip()) * u.kilocalorie_per_mole for line in file]
    #print("tabulated data type = ", type(tabulated_data))
    #print(tabulated_data)
    W_r = mm.Continuous1DFunction(tabulated_data, 1.77500000 * u.angstroms, 50.000000 * u.angstroms)
    MP_potential.addTabulatedFunction('W', W_r)


for i in range(nom):
        MP_potential.addParticle([TypeMP[i],q_bead[i] * Q_DH])
for i in range(nom,ntot):
        MP_potential.addParticle([TypeMP[i],q_bead[i] * u.elementary_charge])

for i in range(0,int(nom/5)):
    iP1 = 5*i    # P1
    iP2 = 5*i+1  # P2
    iS1 = 5*i+2  # S1
    iB  = 5*i+3  # B
    iS2 = 5*i+4  # S2
    MP_potential.addExclusion(iP1,iP2)
    MP_potential.addExclusion(iP2,iS1)
    MP_potential.addExclusion(iS1,iB)
    MP_potential.addExclusion(iS1,iS2)
    if(i < int(nom/5)-1 ):
        MP_potential.addExclusion(iS2,iS2+1)

system.addForce(MP_potential)


#~~~~~~~~~~~~~~~~  Add Bonded P-P Debye-Huckel Interaction using Custom Bonded Interaction ~~~~~~~~~~~~~#
'''
Bond_DHforce = mm.CustomBondForce("BscaleDH * Bq1q2 * exp(-BkappaDH*r)/r")
Bond_DHforce.addGlobalParameter("BscaleDH", l_Bjerrum * u.kilocalorie_per_mole / u.elementary_charge**2)
Bond_DHforce.addGlobalParameter("BkappaDH", kappaDH)
Bond_DHforce.addPerBondParameter('Bq1q2')
for i in range(0,int(nom/5)):
    iP1 = 5*i    # P1
    iP2 = 5*i+1  # P2
    tempBq1q2 = q_bead[iP1]*q_bead[iP2]* Q_DH**2
    Bond_DHforce.addBond(iP1,iP2,[tempBq1q2]) 
system.addForce(Bond_DHforce)
'''

#~~~~~~~~~~~~~~~~  Nonbonded Coulombic Interaction with PME ~~~~~~~~~~~~~#
##^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ EXCLUDED VOLUME INTERACTION BETWEEN PROTEIN-ION AND ION-ION (did not include protein-protein)^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^##
D0 = [
  #    A        U        C        G       S        P=P3     Mg2+    Ca2+   
    0.32000, 0.32000, 0.32000, 0.32000, 0.51300, 0.44100, 0.45200, 0.53200,   # A  
    0.32000, 0.32000, 0.32000, 0.32000, 0.50400, 0.43200, 0.44200, 0.52300,   # U   
    0.32000, 0.32000, 0.32000, 0.32000, 0.50400, 0.43200, 0.44300, 0.52300,   # C
    0.32000, 0.32000, 0.32000, 0.32000, 0.53100, 0.45900, 0.47000, 0.55000,   # G
    0.51300, 0.50400, 0.50400, 0.53100, 0.52200, 0.45000, 0.46100, 0.54100,   # S
    0.44100, 0.43200, 0.43200, 0.45900, 0.45000, 0.37800, 0.38900, 0.46900,   # P=P3
    0.45200, 0.44200, 0.44300, 0.47000, 0.46100, 0.38900, 0.40000, 0.48000,   # Mg2+
    0.53200, 0.52300, 0.52300, 0.55000, 0.54100, 0.46900, 0.48000, 0.56000]   # Ca2+
beadTypeEV = []
for i in range(0,nom):
    if ((i%5) == 0 or (i%5) == 1): #P
            beadTypeEV.append(5)
    elif ((i%5) == 2 or (i%5) == 4): #S
            beadTypeEV.append(4)
    else:
        if (inputfile['Base'] == "A"):
            beadTypeEV.append(0) #B
        elif (inputfile['Base'] == "U"):
            beadTypeEV.append(1) #B
        elif (inputfile['Base'] == "C"):
            beadTypeEV.append(2) #B
        elif (inputfile['Base'] == "G"):
            beadTypeEV.append(3) #B
for i in range(nom+nNeg+nPos,ntot):
    if (inputfile['Ion'] == "Mg"):
        beadTypeEV.append(6)
    elif (inputfile['Ion'] == "Ca"):
        beadTypeEV.append(7)
D0 = u.Quantity(D0,u.nanometers) ## 21 X 21 matrix initialize


switchEV = [
  # P               B/S             Mg 
    1,              1,              0,      # P  
    1,              1,              1,      # B/S
    0,              1,              1]      # Mg2+

TypeEV = []
for i in range(0,nom):
    if ((i%5) == 0 or (i%5) == 1): #P
            TypeEV.append(0)
    else:
            TypeEV.append(1) #B
for i in range(nom,ntot):
    TypeEV.append(2)



energy_function = "switchEV(typeEV1,typeEV2) * step(sig - r) * epsEV * ((R6 - 2) * R6 + 1); R6 = (sig / r)^6; sig = D0(type1, type2);"
EVProtIon_interaction = mm.CustomNonbondedForce(energy_function)
EVProtIon_interaction.addTabulatedFunction('switchEV', mm.Discrete2DFunction(3, 3, switchEV))
EVProtIon_interaction.addPerParticleParameter('typeEV')
EVProtIon_interaction.addTabulatedFunction('D0', mm.Discrete2DFunction(8, 8, D0))  # D0 from Hung PNAS 2019
EVProtIon_interaction.addGlobalParameter('epsEV', 1.0 * u.kilocalories_per_mole)
EVProtIon_interaction.addPerParticleParameter('type')
EVProtIon_interaction.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
EVProtIon_interaction.setCutoffDistance(10.0*unit_len)

for i in range(ntot):
    EVProtIon_interaction.addParticle([TypeEV[i],beadTypeEV[i]])

'''
EVProtIon_interaction = mm.CustomNonbondedForce("epsEVPI * step(sig-r) * (((sig/r)^12-2*(sig/r)^6)+1);epsEVPI=sqrt(epsEVPI1*epsEVPI2);sig=sig1+sig2")
EVProtIon_interaction.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
EVProtIon_interaction.addPerParticleParameter('epsEVPI')
EVProtIon_interaction.addPerParticleParameter('sig')
EVProtIon_interaction.setCutoffDistance(10.0*unit_len)
EVProtIon_interaction.setUseSwitchingFunction(False)

for i in range(ntot):
    tempepsEVPI  = epsEVPI[i] *  mm.KJPerKcal * u.kilojoules_per_mole
    tempsigma    = sigma[i] * unit_len  
    EVProtIon_interaction.addParticle([tempepsEVPI,tempsigma])
    #EVProtIon_interaction.addParticle([beadTypeEV[i], tempepsEVPI])
'''
for i in range(0,int(nom/5)):
        iP1 = 5*i    # P1
        iP2 = 5*i+1  # P2
        iS1 = 5*i+2  # S1
        iB  = 5*i+3  # B
        iS2 = 5*i+4  # S2
        EVProtIon_interaction.addExclusion(iP1,iP2)  # P1-P2
        EVProtIon_interaction.addExclusion(iP2,iS1)  # P2-S1
        EVProtIon_interaction.addExclusion(iS1,iB) # A-S1
        EVProtIon_interaction.addExclusion(iS1,iS2) # S1-S2
        if(i < int(nom/5)-1 ):
                EVProtIon_interaction.addExclusion(iS2,iS2+1) # S2-(P1+1)

system.addForce(EVProtIon_interaction)


    
cm = mm.CMMotionRemover()
system.addForce(cm)    
    
    
for i in range(system.getNumForces()):
    force = system.getForce(i)
    force.setForceGroup(i)


# In[5]:



platform = mm.Platform.getPlatformByName('CUDA')
properties = {'Precision': 'mixed'}
integratorLD = mm.LangevinIntegrator(T,fric_coef,timeStep)
#simulation = app.Simulation(pdb.topology,system,integratorLD,platform)
simulation = app.Simulation(pdb.topology,system,integratorLD,platform,properties)

if inputfile['rstsim'] == "False":
        simulation.context.setPositions(positions)
        simulation.context.setVelocitiesToTemperature(T)
        print("restart: False\n")
else:
        print("restart: True\n")
        simulation.loadCheckpoint('chkin.chk')

simulation.reporters.append(app.DCDReporter('animation.dcd',timeint))
simulation.reporters.append(app.StateDataReporter('property.dat',timeint,potentialEnergy=True,totalEnergy=True,step=True,time=True,temperature=True,speed=True,separator='\t'))
app.PDBFile.writeFile(pdb.topology,positions,open('topology.pdb', 'w'))
simulation.reporters.append(app.CheckpointReporter('chkout.chk', timeint))

print("\nEnergy components of initial position\n")
for i in range(system.getNumForces()):
    print(i,simulation.context.getState(getEnergy=True,groups={i}).getPotentialEnergy())

simulation.step(totstep)

print("\nEnergy components of final position\n")
for i in range(system.getNumForces()):
    print(i,simulation.context.getState(getEnergy=True,groups={i}).getPotentialEnergy())

state = simulation.context.getState(getPositions=True)
finalpos = state.getPositions(asNumpy=True).value_in_unit(u.nanometers)



finalpositions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(pdb.topology, finalpositions, open('finalpos.pdb', 'w'))

time_en = datetime.datetime.now()
print("Simulation stopped at\t")
print(time_en)
print("time taken\t")
print(time_en-time_st)




