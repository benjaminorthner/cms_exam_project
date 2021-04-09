from config import *

####################
# Helper Functions #
####################

# returns only up/down spin part of state vector
def up(x): return x[0:NSITES]
def down(x): return x[NSITES:2*NSITES]

# returns total occpuation of a site (QD site has site index 0)
def n_site(x,site):
    if up(x)[site] == 1 and down(x)[site] == 1: # double occupied
        return 2
    
    elif up(x)[site] == 1 or down(x)[site] == 1: # single occupied
        return 1
    
    return 0 # not occupied

# returns list of binary digits of an integer i, padded on the left with zeroes
def toBinary(i, length):
    bin_i = bin(i)[2:]
    padding = length - len(bin_i)
    
    return np.array(padding*[0] + [int(x) for x in bin_i])


###################
# Generate States #
###################

# Generate all possible states of the system,
# and save states with NUP and NDOWN to a list (occpuation number basis vectors)
def generateStates(NUP, NDOWN):
    states = []
    for i in range(0, 4**NSITES):
        state = toBinary(i, 2*NSITES)
        
        if sum(up(state)) == NUP and sum(down(state)) == NDOWN:
            states += [state]
    
    return states


########################
# Generate Hamiltonian #
########################

# generate pairs of spin sites with allowed hoppings (only between QD and bath with same spin)
hoppings = np.array([[i, j+i] for i in [0, NSITES] for j in range(1,NBATH+1)], dtype=int)
hoppings = np.append(hoppings, np.flip(hoppings,1), 0) # add reverse hoppings
hoppings = hoppings.tolist() # problem later if using numpy array

# function that checks if states have a single hopping between them by
# calculating difference between states
# if true returns the hopping indices
# if false returns -1
# ordering of returned indices is such that the first index represents an
# occupied state in state 1 and the second an occupied state in state2
def hop_check(state1, state2):

        state_diff = state2 - state1

        # checks that only 1 hopping happened
        if sum(abs(state_diff)) == 2:
            # extracts indices of hopping pair
            test_hop = np.array([], dtype=int)
            for i, spinsite in enumerate(state_diff):
                if abs(spinsite) == 1:
                    test_hop = np.append(test_hop, [i*spinsite]) # *spinsite to transfer sign of difference
                test_hop = np.sort(test_hop)                     #  element to index for sorting
            return abs(test_hop).tolist()
        return [-1,-1]


# returns sign of the hopping based on commutation of creation/annihilation operators
# important: order of state1/2 must be the same as it was when computing test_hop
def hop_sign(test_hop, state1, state2):
    # checks no of single commutations needed before annihilation with vacuum and raises -1 to that power
    # harmitian conjugate term always vanishes and thus is not accounted for here. 
    # (Ordering of hopping is designed such that this never vanishes if a hopping exists)
    return (-1)**sum(state1[0:test_hop[0]])*(-1)**sum(state2[0:test_hop[1]])


def generateHamiltonian(states, Vg):

    ##############################
    # Diagonal Hamiltonian Terms #
    ##############################

    # initialize
    dim = len(states)
    H_coulomb = np.zeros((dim,dim))
    H_qd = np.zeros((dim,dim))
    H_bath = np.zeros((dim,dim))

    # generate
    for i,state in enumerate(states):

        # if QD singly occupied
        if n_site(state, 0) == 1:
            H_qd[i][i] += E - e*Vg

        # if QD doubly occupied
        if n_site(state, 0) == 2:
            H_coulomb[i][i] += U
            H_qd[i][i] += 2*(E - e*Vg)

        # if Bath sites occupied
        for bsite in range(1, NBATH+1):

            if n_site(state,bsite) == 2:
                H_bath[i][i] += 2 * Ek[bsite - 1]

            elif n_site(state,bsite) == 1:
                H_bath[i][i] += Ek[bsite - 1]

    ##################################
    # Off-Diagonal Hamiltonian Terms #
    ##################################

    # initialize
    H_bath_lead = np.zeros((dim,dim))

    # generate
    for i, state1 in enumerate(states):
        for j, state2 in enumerate(states):
            test_hop = hop_check(state1,state2)

            if test_hop in hoppings:
                k = np.sort(test_hop)[1]%NSITES - 1 # value of k
                H_bath_lead[i,j] += hop_sign(test_hop, state1, state2) * Vk[k]
                    
    
    return H_qd + H_coulomb + H_bath + H_bath_lead












    ################################################
############### NEW STUFF ############################
    ################################################


def gen_occupation_vector(state, states):
    """
    generates a vector of occupation numbers from a state v
    :param state: np.array in state basis (len(v) == len(states))
    :param states: list of all possible (np.array) states given in occupation number basis 
    :return: vector v_i_s = <n_i_s>
    """

    # initialise vector of occupation numbers
    occ_vector = np.zeros([len(states[0])])

    # we square the weight because one contribution from bra and one from ket
    for i, istate in enumerate(states):
        occ_vector += (abs(state[i])**2)*istate

    return occ_vector


def timeEvolve(state, Hamiltonian, t):
    """
    Returns a state time evolved by t units with a time-independent Hamiltonian
    """
    state_t = expm(1j * Hamiltonian * t).dot(state)
    return state_t


def visualise_occupation(occ_vector, show=False):
    data = np.row_stack([occ_vector[:NSITES], occ_vector[NSITES:]])  
    X = np.arange(NSITES)

    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])

    width = 0.4
    labels = ["QD", "B1", "B2", "B3"]
    ax.set_xticks(X)
    ax.set_xticklabels(labels)
    

    ax.bar(X - width/2, data[0], color = 'b', width = width)
    ax.bar(X + width/2, data[1], color = 'g', width = width)

    if show:
        plt.show()