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
x2 = x1+(10*10**(-6))*(wSI/sc.c) #ca.50+205

#die Verteilung der heißen e- als fkt

def n(x,r):
	k = 2
	e = 2.718281828459045
	if x<x1:
		n=0
	elif x>x2:
		n=0
	else:
		n = 0.5*ne*(e**(-((r+k*0.5*sigma)**2)/(2*sigma**2))+e**(-((r-k*0.5*sigma)**2)/(2*sigma**2)))
	return n

def schicht(x,r):
	e = 2.718281828459045
	x3 = x2 + 2
	if x<x2:
		n=0
	elif x>x3:
		n=0
	else:
		if r<5.*sigma:
			n = 2.5*ne
		else:
			n = 0
	return n

"""
to counteract numerical heating: dx <= Ld
CFL: dt < dx 
"""


#Debyelänge muss aufgelöst werden gegen Gridheating, lambda_d = 11,88
dx = 2. #dx muss kleiner sein damit ich die Protonenschicht auch dünn machen kann
dr = 2. #ist kleiner gegenüber sigma=100	
Lx = 8.*dx*round(750./(8*dx)) #von den bisherigen Teilchen Trajektorien gneommen
Lr =  8.*dr*round(7.5*sigma/(8*dr))            #simulation length
tsim = 630 #630+70              #duration of the simulation
#Höchste verwendete Mode
m = 1
dt = 1/np.sqrt((1+m**2)/dr**2 + 1/dx**2) 
nppcion = 81 #ist 30²
nppc= 441 #8²       
 
diagEveryS = int(tsim/100.)  #frequency of outpu
diagEveryT = int(tsim/30.)
diagEveryF = int(tsim/15.)

# --------------------------------------
# SMILEI's VARIABLES (DEFINED IN BLOCKS)
# --------------------------------------
 
Main(
    geometry = "AMcylindrical",
     
    number_of_AM = m+1,
    
    solve_poisson = False,
     
    timestep = dt,
    simulation_time = tsim,
     
    cell_length = [dx, dr],
    grid_length  = [Lx, Lr],
     
    number_of_patches = [ 8,8 ],#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
     
    EM_boundary_conditions = [['silver-muller', 'silver-muller'],['buneman', 'buneman']] , 

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
    name = 'eonh',
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
    number_density = schicht,
#    is_test = True,
    boundary_conditions = [
        ['remove', 'remove'],
	['remove', 'remove']
    ]#,
#    time_frozen = 70
)


#neutralisierende Elektronen
Species(
    name = 'eon',
    position_initialization = 'random',
    momentum_initialization = 'maxwell-juettner',
    particles_per_cell = nppc,
    mass = 1.,
    charge = -1.0,
    number_density = schicht,
    temperature = [0.08*Te],
#    is_test = True,
    boundary_conditions = [
        ['remove', 'remove'],
	['remove', 'remove']
    ]#,
#    time_frozen = 70
)

LoadBalancing(
    every = 30
)
 
### DIAGNOSTICS
 

#nur Sachen innerhalb einer 5sigma umgebung um die Mitte betrachten alles andere frist unnötig ressourcen
def my_filter(particles):
    return ((particles.y**2+particles.x**2)**0.5<510)*(particles.px>0)

#insgesammt 499 Zeitschritte ich bin nur am Kram am Ende interessiert
DiagTrackParticles(
    species = "pro",
    every = [444,20],#mich interessiert primär der Zustand am Ende
#    flush_every = 100,
#    filter = my_filter,
     attributes = ["x","y","z", "px", "py","pz"]
)

"""
DiagFields(
    every = 80,
    fields = ["El","Er","Et"],
)


DiagFields(
    every = 50,
    fields = ["Rho_eonn","Rho_eonh","Rho_pro"],
)
"""
