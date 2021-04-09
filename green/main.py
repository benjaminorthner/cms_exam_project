from functions import *

##############################################################
## Finding No. of electrons with lowest energy ground state ##
##############################################################

Vg = 0
ground_states = []
for NUP in range(NSITES+1):
    for NDOWN in range(NUP,NSITES+1):
        if NUP == NDOWN == 0: continue

        states = generateStates(NUP, NDOWN)
        Hamiltonian = generateHamiltonian(states, Vg)

        ground_states += [{"Emin":min(np.linalg.eigvals(Hamiltonian)), "up":NUP, "down":NDOWN}]

############################################################
## Find NUP & NDOWN that gives lowest ground state energy ##
############################################################

energies = [gs["Emin"] for gs in ground_states]
imin = energies.index(min(energies))

NUP = ground_states[imin]["up"]
NDOWN = ground_states[imin]["down"]

print("\nLowest ground state energy (Vg=0) for NUP = {} and NDOWN = {}".format(NUP, NDOWN))

#####################################
## (n0 vs Vg) and (g vs Vg) graphs ##
#####################################

states = generateStates(NUP, NDOWN)

#np.random.seed(1)
#randState = np.random.rand(len(groundState))
#randState = NSITES * randState / sum(randState)

Hamiltonian = generateHamiltonian(states, 0)
# extract GS energy and its index in the list
Emin, idx = min((val, idx) for (idx, val) in enumerate(np.linalg.eigvals(Hamiltonian)))
# extract GS
groundState = np.linalg.eig(Hamiltonian)[1][:,idx]

ov = gen_occupation_vector(groundState, states)
ov = gen_occupation_vector(timeEvolve(groundState, Hamiltonian, 10), states)
visualise_occupation(ov, show=True)


