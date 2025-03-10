"""
The State module for use in AtomArray. Store states as States compatible with Experiment and 
initialize particular state in Experiment. 

Author: Emil Henningsen, tzs820, s240670

"""

from utils import *
#from StateCollection import StateCollection #for use in matmul with StateCollection .. unavailable for now

class State:

    #Internal storage:

    system_size = None

    name = None

    value = None
    right_vector = None
    left_vector = None
    hermiticity_metric_single = None    #The W_n from the chinese NH article

    _evolved_right_vector = None
    _evolved_left_vector = None
    _evolved_hermiticity_metric_single = None

    #Private variables:

    #Warnings:
    """Operator overload"""
    def __init__(self, N : int, name : int, val : dtypes, rvec : ndarray, lvec : ndarray = None, compute_left : bool = True):
        """
        Initialize State object. 

        --- INPUT

        N: int - number of atoms (i.e. system size)
        name: int - name of state (usually just the number). Will be used in dictionary storage in Atom Array to keep track of state through different instances of Hamiltonian,
        e.g. when sweeping through a parameter.
        val: dtypes - value of state. value w.r.t. current instance of Hamiltonian. To be updated throughout the instances. See utils for allowed dtypes.
        rvec: ndarray of dtypes of shape (N,) - Right Vector w.r.t current instance of Hamiltonian. Also to be updated.
        lvec (optional): ndarray of dtypes of shape (N,) - left Vector ... Biorthogonal with rvec, i.e. bra(lvec_i)ket(rvec_j) = delta(i,j). Can be obtained from rvec by hermiticity metric: bra(rvec)@P.
        compute_Wn (optional): bool. Toggle whether or not to compute hermiticity metric single, W_n from chinese NH article. Also computes the left vector. Default: True         
        """

        self.system_size = N
        self.name = name
        self.value = val
        self.right_vector = rvec
        self._evolved_right_vector = deepcopy(rvec)

        self.assertDataType(val)
        self.assertVectorProps(rvec)

        if lvec is not None:
            self.assertVectorProps(lvec)
            self.left_vector = lvec
            self._evolved_left_vector = deepcopy(lvec)
        elif compute_left:
            self.computeHermiticityMetricSingle()
            self.computeEvolvedHermiticityMetricSingle()
            self.computeLeftVector()

    #We don't define get and set item methods for this class, because we want to do checks on rvec / lvec, whenever they are changed.

    @singledispatchmethod
    def __matmul__(self, B : Self) -> dtypes:
        """matrix multiplication of states, bra(L_A) @ ket(R_B)"""
        self.assertLeftVectorInit()
        B.assertRightVectorInit()
        return self.left_vector @ B.getRightVector()
    
    @__matmul__.register
    def _(self, B : ndarray) -> ndarray:
        """bra(L_A) @ U = (bra(L_A)*ket(U[:,0]), bra(L_A)*ket(U[:,1]), ...)"""
        self.assertLeftVectorInit()
        #assert(len(self.left_vector) == len(B)) #is this necessary tho? Doesn't @ already do this check? 
        return self.left_vector @ B
    
    @__matmul__.register
    def _(self, B : list) -> list[dtypes]:
        """bra(L_A) @ (ket(R_B_0), ket(R_B_1), ...) braket of list of kets (States)"""
        self.assertLeftVectorInit()
        result = []
        for state in B:
            state.assertRightVectorInit()
            result.append(self.left_vector @ state.getRightVector())
        return result
    
    #@__matmul__.register #TODO: Figure out how to type hint StateCollection without importing... One way could be to do method in StateCollection??? 
    #def _(self, B : StateCollection) -> list[dtypes]:
    #    """Equivalendt to A @ list[B], but with StateCollection instead of list of states"""
    #    #self.assertLeftVectorInit()
    #    list_of_states = list(B.getStates().values())
    #    return self @ list_of_states
    
    """Assertions"""
    def assertDataType(self, x : dtypes) -> None:
        """
        Assert if x is of allowed datatypes. See utils for allowed datatypes.
        """
        assert(x.dtype in get_args(dtypes))

    def assertRightVectorInit(self) -> None:
        assert(self.right_vector is not None)
    
    def assertLeftVectorInit(self) -> None:
        assert(self.left_vector is not None)

    def assertVectorLength(self, vec : ndarray) -> None:
        assert(vec.shape == (self.system_size,))

    def assertVectorNorm(self, vec : ndarray) -> None:
        """
        Assert vector norm less than 1 (+ tolerance, see utils for tolerance level).
        """
        probability = np.sum(np.abs(vec)**2)
        assert(probability < 1+tolerance)

    def assertVectorProps(self, vec: ndarray) -> None:
        self.assertDataType(vec)
        self.assertVectorLength(vec)
        self.assertVectorNorm(vec)

    def assertHermiticityMetricSingleInit(self) -> None:
        assert(self.hermiticity_metric_single is not None)

    """GET methods"""

    def getName(self) -> str:
        return self.name

    def getValue(self) -> dtypes:
        return self.value

    def getDecayRate(self) -> dtypes:
        return -2 * imag(self.value)

    def getEnergyShift(self) -> dtypes:
        return real(self.value)

    def getVector(self) -> ndarray:
        """Default Vector is right Vector"""
        return self.getRightVector()

    def getEvolvedVector(self) -> ndarray:
        return self.getEvolvedRightVector()

    def getRightVector(self) -> ndarray:
        return self.right_vector
    
    def getLeftVector(self) -> ndarray:
        return self.left_vector

    def getEvolvedRightVector(self) -> ndarray:
        return self._evolved_right_vector

    def getEvolvedLeftVector(self) -> ndarray:
        return self._evolved_left_vector

    def getHermiticityMetricSingle(self) -> ndarray:
        """Get the W_n for this state, i.e. inverse(rvec*rvec^{dagger})"""
        return self.hermiticity_metric_single

    """SET methods"""

    def setName(self, name : str) -> None:
        self.name = name

    def setValue(self, val : dtypes) -> None:
        
        #Assertations
        self.assertDataType(val)
        
        self.value = val

    def setVector(self, vec : ndarray) -> None:
        self.setRightVector(vec)

    def setEvolvedVector(self, ev_vec : ndarray) -> None:
        self.setEvolvedRightVector(ev_vec)

    def setRightVector(self, vec : ndarray, compute_left : bool = True) -> None:

        #Assert
        #self.assertVectorProps(vec)
        self.right_vector = vec

        if compute_left:
            self.computeHermiticityMetricSingle()
            self.computeLeftVector()

    def setLeftVector(self, vec : ndarray) -> None:

        #Assert
        #self.assertVectorProps(vec)

        self.left_vector = vec

    def setEvolvedRightVector(self, ev_rvec : ndarray, compute_left : bool = True) -> None:
        self._evolved_right_vector = ev_rvec
        if compute_left:
            self.computeEvolvedHermiticityMetricSingle()
            self.computeEvolvedLeftVector()

    def setEvolvedLeftVector(self, ev_lvec : ndarray) -> None:
        self._evolved_left_vector = ev_lvec

    def setHermiticityMetricSingle(self, Wn : ndarray) -> None:
        self.hermiticity_metric_single = Wn

    """Functions: """

    def signFactor(self, sign : int, evolved : bool = True) -> None:
        self.right_vector *= sign
        self.left_vector *= sign
        if evolved:
            self._evolved_right_vector *= sign
            #self._evolved_left_vector *= sign #TODO: Fix

    def computeHermiticityMetricSingle(self) -> None:
        """
        Compute hermiticity metric, W_n for this state. 

        definition: W_n = [ket(R_n) bra(R_n)]^-1
        """
        self.assertRightVectorInit()
        self.hermiticity_metric_single = outer(self.right_vector, dagger(self.right_vector))    #because vectors, we need to use outer instead of @. 

    def computeEvolvedHermiticityMetricSingle(self) -> None:
        self._evolved_hermiticity_metric_single = outer(self._evolved_right_vector, dagger(self._evolved_right_vector))

    def computeLeftVector(self, rvec : ndarray = None) -> None:
        """
        When initialized, compute the lvec and store. From hermiticity metric. 
        """
        self.assertRightVectorInit()
        self.assertHermiticityMetricSingleInit()
        if rvec is None:
            rvec = self.right_vector
        
        self.left_vector = dagger(rvec) @ self.hermiticity_metric_single

    def computeEvolvedLeftVector(self, ev_rvec : ndarray = None) -> None:
        if ev_rvec is None:
            ev_rvec = self._evolved_right_vector
        #self._evolved_left_vector = dagger(ev_rvec) @ self.hermiticity_metric_single
        self._evolved_left_vector = dagger(ev_rvec) @ self._evolved_hermiticity_metric_single

    def evolveState(self, integrated_eigval : dtypes, evolve_decay : bool = True) -> None:
        """
        evolve eigenvector with integrated decay rate using closed system QM propagator
        """
        if evolve_decay:
            self._evolved_right_vector *= exp(-1j * integrated_eigval)
        else:
            self._evolved_right_vector *= exp(-1j * real(integrated_eigval))

        self.computeEvolvedHermiticityMetricSingle()
        self.computeEvolvedLeftVector()

