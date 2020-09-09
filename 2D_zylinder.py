# -------------------
# MY PYTHON VARIABLES
# -------------------
 
import math 
import numpy as np
import scipy.constants as sc

"""
Physikalischer Aufbau wie bei 10_22_div mitteln
diesmal jedoch ohne Merging mit heißer Debye länge und vergelich von den Fielddiagnostic mit und ohne
Mitteln
"""

ne = 0.083  # 

"""
Parameter aus Lecz S.75
"""

Te = 11.716 #m_e*c² ist Referenztemp, hier die Temperatur für Phelixparameter

"""
Die Referenzzeit (und damit auch die Referenzlänge, wird um den Faktor 
100 verkleinert, da die Intervalle für die Diagnostics >=1 sein müssen.
Müssen die Referenzeinheiten so gewählt werden, dass die Abmessungen der
Simulation in den smilei-Einheiten nicht zu klein sind
taumax und Lx müssen an die neuen Einheiten angepasst werden
"""
#Parameter von lecz S.75
wSI =100* 1.23*10**14 #1/s, ist Plasma Ionenfrequenz da T_r = 1/w_r ist mal 100 eine Verkleinerung von T_r um Faktor 100
sigmaSI = 2.5*10**(-6)#m
sigma = sigmaSI*(wSI/sc.c)#70 L_r

x1 = 50
x2 = x1+(10*10**(-6))*(wSI/sc.c) #ca. 50+410

#die Verteilung der heißen e- als fkt
#die Verteilung der heißen e- als fkt
def n(x,y):
	if x<x1:
		n=0
	elif x>x2:
		n=0
	else:
		if y < -sigma*(3.141592653589793/2)**0.5 + 11*sigma:
			n=0
		elif y > sigma*(3.141592653589793/2)**0.5 + 11*sigma:
			n=0
		else:
			n=ne
	return n

"""
to counteract numerical heating: dx <= Ld
CFL: dt < dx 
"""


#Debyelänge muss aufgelöst werden gegen Gridheating, lambda_d = 11,88
dx = 2. #dx muss kleiner sein damit ich die Protonenschicht auch dünn machen kann
dy = 2. #ist kleiner gegenüber sigma=100	
Lx = 8.*dx*round(750./(8*dx)) #so das auch Teilchen die sich mit c Bewegen nicht abhauen können
Ly =  8.*dy*round(22.*sigma/(8*dy))            #simulation length
tsim = 630         #duration of the simulation
 
nppcion = 9 #ist 30²
nppc= 36 #8²                
 
diagEveryS = int(tsim/100.)  #frequency of outpu
diagEveryT = int(tsim/30.)
diagEveryF = int(tsim/15.)

# --------------------------------------
# SMILEI's VARIABLES (DEFINED IN BLOCKS)
# --------------------------------------
 
Main(
    geometry = "2Dcartesian",
     
    interpolation_order = 2,
     
    timestep_over_CFL = 0.99,
    simulation_time = tsim,
    solve_poisson = False,
     
    cell_length = [dx, dy],
    grid_length  = [Lx, Ly],
     
    number_of_patches = [ 8,8 ],#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
     
    EM_boundary_conditions = [['silver-muller', 'silver-muller'],['silver-muller', 'silver-muller']] , 

    print_every = 10,
     
    random_seed = smilei_mpi_rank
)

Species(
    name = 'ion',
    position_initialization = 'regular',
    momentum_initialization = 'cold',
    particles_per_cell = nppcion,
    mass = 1836., #protonen
    charge = 1.0, # +e ist die Referenzladung
    number_density =  n,
    boundary_conditions = [
        ['remove', 'remove'],
	['remove', 'remove']
    ],
    time_frozen = 2.*tsim
)


Species(
    name = 'eon_h',
    position_initialization = 'regular',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1.0,
    charge = -1.0,
    number_density = n,
    temperature = [Te],
    boundary_conditions = [
        ['remove', 'remove'],
	['remove', 'remove']
    ]
)


#die dicke der Schicht entspricht etwas den 50nm die bei Passoni angegeben werden
#die beschleunigten Protonen als Testteilchen
Species(
    name = 'pro',
    position_initialization = 'random',
    momentum_initialization = 'cold',
    particles_per_cell = nppc,
    mass = 1836.15,
    charge = 1.0,
    number_density = trapezoidal(2.5*ne, xvacuum=x2, xplateau=2, yvacuum=6*sigma, yplateau=10*sigma),
#    is_test = True,	
    boundary_conditions = [
        ['remove', 'remove'],
	['remove', 'remove']
    ]
#    time_frozen = 70
)


#neutralisierende Elektronen
Species(
    name = 'eon_neut',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.0,
    number_density = trapezoidal(2.5*ne, xvacuum=x2, xplateau=2, yvacuum=6*sigma, yplateau=10*sigma),
#    is_test = True,
    temperature = [0.08*Te],
    boundary_conditions = [
        ['remove', 'remove'],
	['remove', 'remove']
    ]
#    time_frozen = 70
)

LoadBalancing(
    every = 40
)
 
### DIAGNOSTICS
 
DiagScalar(
    every = 50
)

"""
DiagParticleBinning(
    #name = "my binning",
    deposited_quantity = "weight",
    every = [550,25],
    time_average = 1,
    species = ["eon_neut", "eon_h"],
    axes = [
        ["px", 0.01, 300, 100],
	["py", 0.01, 150, 50]
    ]
)
"""

DiagTrackParticles(
    species = "pro",
    every = [379,10],#mich interessiert primär der Zustand am Ende
#    flush_every = 100,
#    filter = my_filter,
     attributes = ["x","y", "px", "py","Ex","Ey"]
)

"""
DiagFields(
    every = 80,
    fields = ["Ex","Ey"],
)
"""
"""
Checkpoints(
    # restart_dir = "dump1",
    dump_step = 200,
    dump_minutes = 120.,
    exit_after_dump = False,
    keep_n_dumps = 2,
) 
""" 
