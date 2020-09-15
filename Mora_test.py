# -------------------
# MY PYTHON VARIABLES
# -------------------
 
import math 
import numpy as np
import scipy.constants as sc

ne = 0.186  #initial electron density 

x1 = 0
x2 = 250

"""
Parameter aus Lecz S.75
"""

Te = 380*(10**3)*sc.e/(sc.m_e*sc.c**2) #mc² ist Referenztemp

"""
Die Referenzzeit (und damit auch die Referenzlänge, wird um den Faktor 
100 verkleinert, da die Intervalle für die Diagnostics >=1 sein müüsen.
Müssen die Referenzeinheiten so gewählt werden, dass die Abmessungen der
Simulation in den smilei-Einheiten nicht zu klein sind
taumax und Lx müssen an die neuen Einheiten angepasst werden
"""
wSI =100* 8.28*10**13 #1/s,
sigmaSI = 2.5*10**(-6)#m
sigma = sigmaSI*(wSI/sc.c)

#die Verteilung der heißen e- als fkt
def n(x,y):
	if x<x1:
		n=0
	elif x>x2:
		n=0
	else:
		n = ne*np.exp(-((y-5*sigma)**2)/(2*sigma**2))
	return n

#(hot) debye Länge
lSI = 7.28*10**(-8) #m
l = lSI*(wSI/sc.c)
#gültigkeitsdauer für die isothermen Formeln, ist so direkt in der richtigen Einheit
taumax = 15*np.sqrt(50*(sigmaSI/lSI) + (sigmaSI/lSI)**2)#T_r=807

"""
to counteract numerical heating: dx <= Ld
"""

dx = l		#Debyelänge muss aufgelöst werden gegen Gridheating
Lx = 8.*dx*round(350./(8*dx))
Ly =  8.*dx*round(10.*sigma/(8*dx))            #simulation length
tsim = taumax               #duration of the simulation
 
nppc= 225 #15²                
 
diagEveryS = int(taumax/100.)  #frequency of output
diagEveryF = int(taumax/40.)

# --------------------------------------
# SMILEI's VARIABLES (DEFINED IN BLOCKS)
# --------------------------------------
 
Main(
    geometry = "2Dcartesian",
     
    interpolation_order = 4,
     
    timestep_over_CFL = 0.99,
    simulation_time = tsim,
     
    cell_length = [dx, dx],
    grid_length  = [Lx, Ly],
     
    number_of_patches = [ 8,8 ],
     
    EM_boundary_conditions = [['reflective', 'silver-muller'],['silver-muller', 'silver-muller']] , 

    print_every = 10,
     
    random_seed = smilei_mpi_rank,
    
)

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1836., #protonen
    charge = 1.0, # +e ist die Referenzladung
    number_density =  n,
    boundary_conditions = [
        ['reflective', 'remove'],
	['remove', 'remove']
    ],
#    time_frozen = 2.*tsim
    temperature = [Te]
)
 
Species(
    name = 'eonh',
    position_initialization = 'regular',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1.0,
    charge = -1.0,
    number_density = n,
    temperature = [Te],
    boundary_conditions = [
        ['reflective', 'remove'],
	['remove', 'remove']
    ]
)

 
LoadBalancing(
    every = 20
)
 
### DIAGNOSTICS
 
DiagScalar(
    every = diagEveryS
)
 

DiagFields(
    every = diagEveryS,
    fields = ["Rho_ion","Rho_eonh","Ex"],
    time_average = diagEveryS
)
