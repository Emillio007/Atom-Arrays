"""
Author: tzs820, s240670, Emil Henningsen
Module for constructing different Hamiltonians used for describing atomic arrays interacting with light.
"""

from utils import *
from GreensTensor import *
from Lattice import *

class Hamiltonian:

    #flags:
    initialized = False
    decomposed = False

    #Internal storage:
    N = None
    hamiltonian = None
    eigvec = None
    eigval = None
    index_sorted = None

    #private variables:
    #Manually set diagval, see meeting notes 7/3:
    #diagval = 1 - 1j/2
    #ATTENTION: setting to 0 - 1j/2 instead, as Klaus suggests:
    diagval = - 1j/2

    def __init__(self, N : int = None, ham : ndarray = None):
        if ham == None:
            self.initialized = False
        else:
            self.initialized = True
        self.decomposed = False
        self.hamiltonian = ham
        self.N = N

    #Assertations:
    def assertInit(self) -> None:
        if self.initialized == False:
            raise self.e
        else:
            pass

    def assertDecomposed(self) -> None:
        if self.decomposed == False:
            #push warning
            warnings.warn(warning_hamiltonian_eigendecomp)
            #Then try eigen decomposition:
            self.eigenDecomposition()
        else:
            pass

    """
    GET methods:
    """

    def getN(self) -> int:
        return self.N

    def getHam(self) -> ndarray:
        """ambiguous with AtomArray method, therefore also called getHamMatrix, but kept for backwards compatibility"""
        self.assertInit()
        return self.hamiltonian
    
    def getHamMatrix(self) -> ndarray:
        return self.getHam()

    def getSortedIndex(self) -> ndarray:
        self.assertDecomposed()
        return self.index_sorted

    #Specialized:

    def getEigenDecomp(self, sort : bool = True) -> tuple[ndarray, ndarray]:
        """
        TODO: Description

        sort is not implemented yet.
        """

        #Assertations:
        self.assertInit()
        self.assertDecomposed()

        return self.eigval, self.eigvec

    def getDecayRates(self, sort : bool = True) -> ndarray:
        """
        TODO: Description
        
        Note: If sort = True, the decay rates will be sorted lowest to highest. 
        The sorted indexes are stored if needed for extracting the corresponding eigvec.
        """
        from numpy import imag, argsort, min

        #Reset sorted indexes:
        self.index_sorted = None

        #if ham is empty
        self.assertInit()
        
        #if eigval is empty
        self.assertDecomposed()
        
        decay_rates = -2 * imag(self.eigval)     #Factor of 2, see Asenjo-Garcia et al. NOTE: Any trouble with minus-sign?

        if sort:
            sortind = argsort(decay_rates)
            decay_rates = decay_rates[sortind]
            self.index_sorted = sortind

        """Don't do this: """
        #minval = min(decay_rates)
        #if rescale:
        #    decay_rates -= minval

        return decay_rates
    
    def getAmplNorm(self, index : int) -> ndarray:
        """
        Get amplitude norms for eigenvector at index.

        TODO: Description
        """
        from numpy import sum, imag, real, sqrt, power

        #assertations as getDecayRates:
        self.assertInit()
        self.assertDecomposed()

        if self.index_sorted is None:
            self.getDecayRates(sort=True)

        ind_sort = self.index_sorted

        ampl = self.eigvec[:,ind_sort[index]]
        amplnorm = sqrt(real(ampl)**2 + imag(ampl)**2)  #Get norm of complexvalued 

        accProb = sum(power(amplnorm, 2))
        if accProb > 1:
            warnings.warn(warning_hamiltonian_probabilities)
        print("The probabilities sum to: ", accProb)

        return amplnorm
        

    def isDecomposed(self) -> bool:
        return self.decomposed
    
    def isInit(self) -> bool:
        return self.initialized

    """
    SET methods:
    """

    def setN(self, N : int) -> None:
        self.N = N

    def setHam(self, ham : ndarray) -> None:
        """ambiguous with AtomArray method, therefore also called setHamMatrix, but kept for backwards compatibility"""
        self.hamiltonian = ham

        #flags:
        self.initialized = True
        self.decomposed = False

    def setHamMatrix(self, ham : ndarray) -> None:
        self.setHam(ham)

    """----- Construct standard hamiltonians -----"""

    """
    In this section, the full Hamiltonian is constructed (full 2^N x 2^N space), which is very memory-demanding for larger N. 

    N.B.: Currently incompatible
    """

    def coherence_operators(self, i, j, N):
        """
        Construct the space of coherence operators at i & j and identity elsewhere.
        Ref.: Asenjo-Garcia et al. equation (5)

        --- Parameters

        i: integer, i'th subspace of raising operator
        j: integer, j'th subspace of lowering operator
        N: integer, number of subspaces

        --- Return

        space: Qobj, array (2^N x 2^N) of 0's and 1's representing constructed space.
        """
        from qutip import Qobj, tensor, qeye

        sigma_ge = Qobj([[0, 1], [0, 0]])                         #deexcitation
        sigma_eg = Qobj([[0, 0], [1, 0]])                         #excitation

        if i == 0:
            space = sigma_eg
        elif j == 0:
            space = sigma_ge
        else:
            space = qeye(2)
        for k in range(1, N):
            if k == i:
                space = tensor([space, sigma_eg])
            elif k == j:
                space = tensor([space, sigma_ge])
            else:
                space = tensor([space, qeye(2)])
            #print(k)

        return space

    def H_eff(self, N, w0, D, G):
        """
        Function for calculating effective Hamiltonian, ref.: Asenjo-Garcia et al. equation (5)

        --- Parameters:

        N: integer, number of dipoles (subspaces)
        w0: float, transition frequency of dipoles
        D: array (1 x 3) of complex floats, column-vector of dipole transition elements
        G: array (N x N x 3 x 3) of complex floats, Green's Tensors for a given lattice of dipoles.

        --- Return:

        H_eff: Qobj, array (2^N x 2^N) of complex floats
        """
        from qutip import Qobj
        from scipy.constants import mu_0

        H_eff = 0                                                       #Effective Hamiltonian of the system
        for i in range(N):                                              #Fill the effective Hamiltonian
            for j in range(N):
                if i == j:
                    continue
                else:
                    H_eff += (-mu_0 * w0**2) * D.trans() * Qobj(G[i,j]) * D * self.coherence_operators(i, j, N)
        return H_eff

    def H(self, N, w0, H_eff):
        """
        Function for calculating full Hamiltonian, ref.: Asenjo-Garcia et al. equation (5)

        --- Parameters:

        N: integer, number of dipoles (subspaces)
        w0: float, transition frequency of dipoles
        H_eff: array (2^N x 2^N) of complex floats, effective Hamiltonian

        --- Return:

        H: Qobj, array (2^N x 2^N) of complex floats, full Hamiltonian.
        """
        from qutip import Qobj, qeye, tensor
        from scipy.constants import hbar

        sigma_ee = Qobj([[0, 0], [0, 1]])                             #excited state
        H = 0
        for i in range(N):

            if i == 0:
                space = sigma_ee
            else:
                space = qeye(2)
            for k in range(1, N):
                if k == i:
                    space = tensor([space, sigma_ee])
                else:
                    space = tensor([space, qeye(2)])

            H += hbar * w0 * space
            
        H += H_eff
        return H

    """
    BLOCK HAMILTONIANS IN BASIS OF {|e_j>} - denoting the single-excitation states of excitation at j'th atom and ground-state elsewhere.
    The Hamiltonian as described in e.g. Asenjo-Garcia commute with the number-operator meaning it conserves excitation-numbers. 
    In other words, the Hamiltonian will be of block-structure nicely divided in 0, 1, 2, 3, ... and so on excitations. 

    Result is arrays of size e.g. (N x N) in the single-excitation case. 
    """

    def block(self, N : int, G : ndarray, n : ndarray) -> None:
        """
        Desciption: TODO
        """
        from numpy import zeros
        from scipy.constants import pi

        #Raise flags:
        self.initialized = False
        self.decomposed = False

        self.N = N

        #If only one polarization vector is passed, e.g. ez, then fill (N x 3) dim array of given polarization
        if n.ndim == 1:
            arr = zeros((N,3))
            for i in range(N):
                arr[i] = n
            n = arr

        #We want Hamiltonian in basis of {|e_j>}, meaning one excitation on j'th subspace.
        #Do this by computing every matrix element of NxN matrix. It is the block of one excitation in full Hamiltonian
        block = zeros((N,N), dtype=complex)

        for i in range(N):
            for j in range(N):
                if i == j:
                    block[i, j] = self.diagval
                else:
                    block[i, j] += -3*pi * n[i].transpose() @ G[i,j] @ n[j]        #Dipoles polarized along z-direction.  
                    """IMPORTANT NOTE: 
                    Has been changed to n[i] @ G[i,j] @ n[j], as this is actually the correct equation. 
                    I.e. the field at i'th dipole resulting from j'th dipole. However, it should have no effect
                    on the cases, where the dipoles are all polarized in the same direction. BUT, it might 
                    have an effect on the cases, where they are oriented individually... Struggles.
                    """

        #Store in internal container.
        self.hamiltonian = block #finish with a minus sign due to eigenvalues being damping. See meeting notes 27/3
        #NOTE: 14/10, removed the minus sign, as it turns out to be correct... ? what the hell happened in spring

        #Flag:
        self.initialized = True

    def blockFromLattice(self, lat : Lattice, G : Union[ndarray, greenstensor_types], w0 : float = 1., dimensionless : bool = True):
        """
        Bloch hamiltonian from lattice object and G
        TODO: Description
        """

        #Extract information from lattice object:
        N = lat.getN()
        rij = lat.getDisplacements()
        polarizations = lat.getPolarizations()

        greens_tensor = None

        #If G is not specified, fill G:
        if type(G) is ndarray:
            greens_tensor = G
        else:
            greens_tensor = fill_G(N, rij, w0, dimensionless, type=G)

        #Do the block Hamiltonian
        self.block(N, greens_tensor, polarizations)

    def scalarham(self, N : int, rij : ndarray, w0 : float = 1, dimensionless : bool = True) -> None:
        """
        See scalar() under GreensTensor.py
        """
        from GreensTensor import scalar
        from scipy.constants import pi

        #Flags:
        self.initialized = False
        self.decomposed = False

        self.setN(N)

        h = scalar(N, rij, w0, dimensionless)

        #Multiply with constants, see notes meetin 7/3
        for i in range(N):
            for j in range(N):
                if i == j:
                    h[i, j] = self.diagval
                else:
                    h[i, j] = -3 * pi * h[i, j]

        #store in internal
        self.hamiltonian = -h   #Finish with a minus sign due to eigenvalues being damping. See meeting notes 27/3
        
        #Lower flag:
        self.initialized = True

    def eigenDecomposition(self) -> None:
        """
        TODO: Description
        """
        from numpy.linalg import eig    #import eig, since ham might not be herm

        #Check if ham is empty
        self.assertInit()
        
        eigval, eigvec = eig(self.hamiltonian)

        #Store:
        self.eigval = eigval
        self.eigvec = eigvec

        #Flags:
        self.decomposed = True
        