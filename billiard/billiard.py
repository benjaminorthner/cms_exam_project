from re import template
import sympy as sym
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import copy
import time
import os


class Billiard:
    '''
    USES THE "EXPANSION METHOD FOR STATIONARY STATES OF QUANTUM BILLIARDS"
    TO COMPUTE EIGENERGIES AND EIGENSTATES OF QUANTUM BILLIARD SYSTEMS

    Implementation of the algorithm described in: https://aapt.scitation.org/doi/10.1119/1.19208

    Note: the stadium shape has not been fully generalised for arbitrary dimensions
    a2 = 2 is required
    '''

    # shape = "quartercircle", "circle", "stadium"
    def __init__(self, shape, a1, a2, M0, V0) -> None:
        self.x1, self.x2 = sym.symbols("x1 x2")
        self.m1, self.m2 = sym.symbols("m1 m2", integer=True, nonzero=True)
        self.n1, self.n2 = sym.symbols("n1, n2", integer=True, nonzero=True)
        self.shape = shape
        self.a1 = a1
        self.a2 = a2
        self.M0 = M0
        self.V0 = V0
        self.A = self.get_area()
        self.L = self.get_circumfrence()

        self.phi_m = sym.sqrt(2.0/self.a1)*sym.sin(sym.pi * self.m1 * self.x1 / self.a1) * sym.sqrt(2.0/self.a2) * sym.sin(sym.pi * self.m2 * self.x2 / self.a2)
        self.phi_n = sym.sqrt(2.0/self.a1)*sym.sin(sym.pi * self.n1 * self.x1 / self.a1) * sym.sqrt(2.0/self.a2) * sym.sin(sym.pi * self.n2 * self.x2 / self.a2)

        if self.shape == 'circle' or self.shape == 'stadium':
            # origin needs to be shifted to center of rectangle
            self.phi_m = self.phi_m.subs([(self.x1, self.x1 - self.a1/2), (self.x2, self.x2 - self.a2/2)])
            self.phi_n = self.phi_n.subs([(self.x1, self.x1 - self.a1/2), (self.x2, self.x2 - self.a2/2)])

        # for calculating the wavefunction later
        self.phi = sym.lambdify([self.x1, self.x2, self.m1, self.m2], self.phi_m)

        self.generate_integer_pairs()

    ###################
    #     V_MATRIX    #
    ###################

    def analytic_integral(self):
        if self.shape == 'quartercircle':
            self.v_analytic = sym.integrate(self.phi_m*self.phi_n, (self.x2, sym.sqrt(1-self.x1**2), 1))
            self.v_analytic = sym.lambdify([self.x1, self.m1, self.m2, self.n1, self.n2], self.v_analytic)

        if self.shape == 'circle':
            v_analytic_top = sym.integrate(self.phi_m*self.phi_n, (self.x2, sym.sqrt(1-self.x1**2), 1))
            self.v_analytic_top = sym.lambdify([self.x1, self.m1, self.m2, self.n1, self.n2], v_analytic_top)

        if self.shape == 'stadium':
            v_analytic_top = sym.integrate(self.phi_m*self.phi_n, (self.x2, sym.sqrt(1-(self.x1 - (self.a1/2 - 1))**2), 1))
            self.v_analytic_top = sym.lambdify([self.x1, self.m1, self.m2, self.n1, self.n2], v_analytic_top)


    # numerical integral
    def vnm(self, m1_, m2_, n1_, n2_):
        if self.shape == 'quartercircle':
            return quad(lambda x: self.v_analytic(x, m1_, m2_, n1_, n2_), 0, 1)[0]
        
        elif self.shape == 'circle':
            if (m1_ + m2_+ n1_ + n2_ )%2 == 0 and (m1_ + n1_)%2 == 0:
                # integral in all four corners is the same (only sign may differ)
                return 4*quad(lambda x: self.v_analytic_top(x, m1_, m2_, n1_, n2_), 0, self.a1/2)[0]
            else:
                return 0

        elif self.shape == 'stadium':
            if (m1_ + m2_+ n1_ + n2_ )%2 == 0 and (m1_ + n1_)%2 == 0:
                # integral in all four corners is the same (only sign may differ)
                return 4*quad(lambda x: self.v_analytic_top(x, m1_, m2_, n1_, n2_), self.a1/2 - 1, self.a1/2)[0]
            else:
                return 0

    def generate_v_matrix(self):
        # check if file already exists:
        foldername = "vmatrices"
        filename = "vmatrix_{}_M{}.txt".format(self.shape, int(self.M0))

        # if exists import vmatrix
        if filename in os.listdir("vmatrices"):
            self.v_matrix = np.loadtxt(foldername + "/" + filename)
            print("V-matrix Imported")

        # else compute v_matrix
        else:
            print("Generating V-matrix")
            start = time.time()

            self.v_matrix = np.zeros([self.M0,self.M0])
            for i,m in enumerate(self.pairs):

                # estimate remaining time
                if i== 0:
                    print("{}/{}\t ETR: tbd".format(i+1, self.M0))
                else:
                    elapsed = time.time() - start
                    time_remaining = sum([self.M0 - x for x in range(i, self.M0)])*elapsed / float(sum([self.M0 - x for x in range(i)]))
                    print("{}/{}\t ETR: {}:{:02d}".format(i+1, self.M0, int(time_remaining // 60), int(time_remaining%60)))

                for j,n in enumerate(self.pairs):
                    
                    if j >= i:
                        self.v_matrix[i,j] = self.vnm(m[0], m[1], n[0], n[1])
                    
                    # because v_matrix is symmetric
                    elif j < i:
                        self.v_matrix[i,j] = self.v_matrix[j, i]
            

            print("Finished generating V-matrix")

            # SAVE V MATRIX TO FILE!
            np.savetxt(foldername + "/" + filename, self.v_matrix)
            


    ################################
    #    HAMILTONIAN & EIGENSTUFF  #
    ################################

    def generate_hamiltonian(self):
        self.hamiltonian = np.zeros([self.M0, self.M0])

        for i,m in enumerate(self.pairs):
            self.hamiltonian[i, i] = np.pi ** 2 * (m[0] ** 2 + (self.a1*m[1]/self.a2) ** 2)

        self.hamiltonian = self.hamiltonian + self.V0 * self.v_matrix

        # calc eigenvals and vecs
        self.eigenvals, self.eigenvecs = np.linalg.eig(self.hamiltonian)

        self.eigenvecs = np.transpose(self.eigenvecs)
        # sort based on ascending eigenvals
        self.eigenvecs = np.array([vec for _, vec in sorted(zip(self.eigenvals, self.eigenvecs))])

        # get eigenvals in different units to match Kaufman, Kosztin, Schulten Paper
        self.kn_evals = np.sort(np.sqrt(self.eigenvals)/self.a1)

        self.En_evals = np.sort(self.kn_evals**2 * self.A/(4*np.pi)) / self.a1**2
    

    ##################################
    #     PLOTTING WAVEFUNCTION      #
    ##################################

    def evalute_wavefunction_grid(self, n, x_pixels, y_pixels):
        if self.shape == 'quartercircle':
            self.X1 = np.linspace(0, 1, x_pixels)
            self.X2 = np.linspace(0, 1, y_pixels)

        elif self.shape == 'circle' or self.shape == 'stadium':
            self.X1 = np.linspace(-self.a1/2, self.a1/2, x_pixels)
            self.X2 = np.linspace(-self.a2/2, self.a2/2, y_pixels)

        self.wavefunction = np.zeros([x_pixels, y_pixels])
        for i, x in enumerate(self.X1):
            for j, y in enumerate(self.X2):
                self.wavefunction[i, j] = abs(self.psi(n, x,y))
                

    def plot_wavefunction(self, n, x_pixels, y_pixels, save=False):

        self.evalute_wavefunction_grid(n, x_pixels, y_pixels)

        padding = 0

        if self.shape == 'quartercircle':
            # plot quartercircle
            plt.plot(np.linspace(0,1,200)*x_pixels, y_pixels*np.sqrt(1-np.linspace(0,1,200)**2), color="black")
            plt.plot([0,0], [0,x_pixels], color="black")
            plt.plot([0,y_pixels], [0,0], color="black")

            # plot wavefuction
            plt.imshow(self.wavefunction, cmap='gray_r', origin="lower", interpolation="bicubic")
            plt.xlim([-padding,x_pixels + padding])
            plt.ylim([-padding,y_pixels + padding])
            plt.axis("off")
        
        elif self.shape == 'circle': 
            # plot circle
            X = np.linspace(0, x_pixels, 200)
            plt.plot(X, y_pixels//2 +x_pixels//2 * np.sqrt(1- ((X-x_pixels//2)/(x_pixels//2))**2), color="black")
            plt.plot(X, y_pixels//2 -x_pixels//2 * np.sqrt(1- ((X-x_pixels//2)/(x_pixels//2))**2), color="black")

            # plot wavefunction
            plt.imshow(self.wavefunction, cmap='gray_r', origin="lower", interpolation="bicubic")
            print(self.wavefunction.shape)
            plt.xlim([-padding,x_pixels + padding])
            plt.ylim([-padding,y_pixels + padding])
            plt.axis("off")

        elif self.shape == 'stadium':
            # plot stadium caps
            Xleft = np.linspace(0, y_pixels/2, 100)
            Xright = np.linspace(3*x_pixels/4, x_pixels, 100)
            # left cap
            plt.plot(Xleft, y_pixels//2 +y_pixels//2 * np.sqrt(1- ((Xleft-y_pixels//2)/(y_pixels//2))**2), color="black")
            plt.plot(Xleft, y_pixels//2 -y_pixels//2 * np.sqrt(1- ((Xleft-y_pixels//2)/(y_pixels//2))**2), color="black")
            # right cap
            plt.plot(Xright, y_pixels//2 +y_pixels//2 * np.sqrt(1- ((Xright-3*x_pixels//4)/(y_pixels//2))**2), color="black")
            plt.plot(Xright, y_pixels//2 -y_pixels//2 * np.sqrt(1- ((Xright-3*x_pixels//4)/(y_pixels//2))**2), color="black")
            # flat parts
            plt.plot([y_pixels/2, 3*x_pixels/4],[y_pixels, y_pixels], color="black", linewidth=1.5)
            plt.plot([y_pixels/2, 3*x_pixels/4],[0, 0], color="black", linewidth=1.5)

            # plot wavefunction
            plt.imshow(self.wavefunction.transpose(), aspect='auto', cmap='gray_r', origin="lower", interpolation="bicubic")

            plt.xlim([-padding,x_pixels + padding])
            plt.ylim([-padding,y_pixels + padding])
            plt.axis("off")

        #plt.title("n = {}".format(n+1))

        if save:
            plt.savefig("plots/WF_{}_{}.png".format(self.shape, n), dpi=100, bbox_inches='tight')

        plt.show()


    def plot_wavefunction3D(self, n, x_pixels, y_pixels):

        self.evalute_wavefunction_grid(n, x_pixels, y_pixels)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        X, Y = np.meshgrid(self.X1, self.X2)
        ax.plot_surface(X,Y, self.wavefunction, cmap=cm.viridis, linewidth=0, antialiased=True)
        fig.show()

        
    ##################################
    #       PLOTTING ENERGIES        #
    ##################################

    def plot_e_energies(self, log=False):
        plt.scatter(range(len(self.En_evals)), self.En_evals, s=6, color="black")
        
        if log:
            plt.yscale("log")
        plt.show()

    ####################
    # HELPER FUNCTIONS #
    ####################

    def compute_everything(self) -> None:
        self.analytic_integral()
        self.generate_v_matrix()
        self.generate_hamiltonian()

    def get_area(self) -> float:
        if self.shape == 'quartercircle':
            return np.pi / 4

        elif self.shape == 'circle':
            return np.pi

        elif self.shape == 'stadium':
            return np.pi + self.a2 * (self.a1 -2)

    def get_circumfrence(self) -> float:
        if self.shape == 'quartercircle':
            return 2 + np.pi/2
        
        elif self.shape == 'circle':
            return 2*np.pi

        elif self.shape == 'stadium':
            return 2*(np.pi + (self.a1-2))

    # generates list of i_max positive integer pairs
    def generate_integer_pairs(self):
        pair = np.array([1,1], dtype=int)

        odd = lambda x: x%2 == 1
        even = lambda x: x%2 == 0

        self.pairs = [copy.deepcopy(pair)]
        
        for _ in range(self.M0-1):
            if odd(pair[0]) and pair[1]==1 :
                pair += np.array([1,0])

            elif pair[0]==1 & even(pair[1]):
                pair += np.array([0,1])

            elif (odd(pair[0]) and even(pair[1])) or (even(pair[0]) and odd(pair[1])):
                pair += np.array([-1,1])

            elif (odd(pair[0]) and odd(pair[1])) or (even(pair[0]) and even(pair[1])):
                pair += np.array([1,-1])
            
            self.pairs.append(copy.deepcopy(pair))
        
    # returns nth wavefuntion at coordinate (x1_, x2_)
    def psi(self, n_, x1_, x2_):

        psi = 0
        # pairs not sorted because order in rows of eigenvecs is according to og arrangement
        for i, m in enumerate(self.pairs):
            psi += self.eigenvecs[n_][i] * self.phi(x1_, x2_, m[0], m[1])
        
        return psi