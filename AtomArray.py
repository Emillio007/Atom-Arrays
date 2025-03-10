"""
The AtomArray class. Used to specify an atom array system with chosen Hamiltonian, vary parameters and evolve in time. 

AUTHOR: Emil Henningsen, tzs820, s240670
"""

#modules:
from Hamiltonian import *
from GreensTensor import *
from Lattice import *
from Parameter import *
from ParameterSet import *
from State import *
from StateCollection import *
from utils import *

class AtomArray:

    """flags:"""
    initialized_lat = False
    initialized_H = False


    """Internal storage:"""
    parameters = None
    num_states = None

    #G = None #maybe this one actually does not belong here?
    lat = None
    H = None

    #hermiticity metric
    P = None
    
    eigen_states : StateCollection = None

    """Private variables:"""

    """INIT: """

    def __init__(self, 
                 parameters : ParameterSet, 
                 lattice : Union[Lattice, lattice_types],
                 hamiltonian : Union[Hamiltonian, hamiltonian_types],
                 ):
        """
        Initialize new AtomArray object. Wraps the functionality of Lattice, GreensTensor and Hamiltonian with Parameter and State. To be used in Experimenter.

        --- INPUT

         - parameters (ParameterSet) - Typed Dictionary of Parameters to describe system. Can be static and/or dynamic. If dynamic, system is initialized in first value of parameters.
         - lattice (Lattice or str) - Lattice or Literal[lattice_types], see utils. The system lattice.
         - hamiltonian (Hamiltonian or str) - Hamiltonian or Literal[hamiltonian_types], see utils. The system Hamiltonian. 
        """
        #Handle parameters

        self.assertParametersStandard(parameters)
        
        self.parameters = parameters

        """
        TODO: Implement way to handle dynamic parameters already now, so it doesn't have to be special case handling in every lattice/hamiltonian update!
        
        """

        #Handle lattice
        #init lattice by updating lattice as auto will get index=0
        self.updateLattice(lattice)
        
        #Handle hamiltonian
        #same as for lattice
        self.updateHamiltonian(hamiltonian)
        
        #Set number of states dependent on parameters:
        self.setNumStates(parameters)

        #initialize eigenstates (right eigvecs) and hermiticity metric:
        self.initEigenStates(sort=True)
            #By default sort lowest decay rate to highest as naming
        
    """ Helper functions """

    def updateLattice(self, lattice : Union[Lattice, lattice_types]) -> None:
        """
        update lattice with current iteration of parameter
        """

        self.initialized_lat = False

        if type(lattice) is Lattice:
            self.lat = lattice
            self.initialized_lat = True
        else:
            #initialize empty lattice object
            self.lat = Lattice()
            try:
                match(lattice):

                    case "linear":
                        if "angle" in self.parameters:
                            if self.parameters["angle"].getStatic():
                                self.lat.linlat(self.parameters["N"]["auto"], 
                                                self.parameters["distance"]["auto"], 
                                                self.parameters["direction"]["auto"], 
                                                self.parameters["polarization"]["auto"])
                            else:
                                angle = self.parameters["angle"]["auto"]
                                self.lat.linlat(self.parameters["N"]["auto"], 
                                                self.parameters["distance"]["auto"], 
                                                self.parameters["direction"]["auto"], 
                                                polarizations=array([cos(angle), 0, sin(angle)]))
                        else:
                            self.lat.linlat(self.parameters["N"]["auto"], 
                                            self.parameters["distance"]["auto"], 
                                            self.parameters["direction"]["auto"], 
                                            self.parameters["polarization"]["auto"])
                        self.initialized_lat = True
                    
                    case "broken":
                        self.lat.twopiece(self.parameters["N"]["auto"], 
                                          self.parameters["distance"]["auto"], 
                                          self.parameters["broken_angle"]["auto"], 
                                          self.parameters["direction"]["auto"], 
                                          self.parameters["polarization"]["auto"])
                        self.initialized_lat = True

                    case "circular":
                        self.lat.circlelat(self.parameters["N"]["auto"], 
                                           self.parameters["distance"]["auto"], 
                                           self.parameters["circle_measure"], 
                                           self.parameters["circle_polarization"], 
                                           self.parameters["polarization"])
                        self.initialized_lat = True

                    case _:
                        print("when trying lattice: ", lattice)
                        raise exception_not_yet_implemented
                    
            except:
                raise exception_required_parameters

    def updateHamiltonian(self, hamiltonian : Union[Hamiltonian, hamiltonian_types]) -> None:
        """
        Update Hamiltonian with current parameter values in lattice
        """

        self.initialized_H = False

        if type(hamiltonian) is Hamiltonian:
            self.H = hamiltonian
            self.initialized_H = True
        else:
            #initialize empty hamiltonian object
            self.H = Hamiltonian()
            try:
                match(hamiltonian):

                    case "block":
                        #TODO: Implement more than 1-excitation manifold functionality. This must also be implemented in block().
                        #TODO: Implement full functionality of w0...
                        self.H.blockFromLattice(self.lat, self.parameters["greenstensor_type"])
                        self.initialized_H = True

                    case _:
                        print("when trying hamiltonian: ", hamiltonian)
                        raise exception_not_yet_implemented
                    
            except:
                raise exception_required_parameters

    """Assertions: """

    def assertParametersStandard(self, parameters) -> None:
        try:
            #Make sure standard parameters are specified by checking type. If excepted, they are not specified correctly, if at all.
            assert(parameters["N"].getDataType() in get_args(dtypes_int))
            assert(parameters["num_excitations"].getDataType() is int or int64 and parameters["num_excitations"][0] == 1)   #For now, only 1-excitation manifold allowed.
            assert(parameters["greenstensor_type"] in get_args(greenstensor_types))
        except:
            raise exception_standard_parameters

    def assertInitializedHamiltonian(self) -> None:
        assert(self.initialized_H)
    
    def assertInitializedLattice(self) -> None:
        assert(self.initialized_lat)

    def assertHamMatrixProperties(self, ham : Hamiltonian = None) -> None:
        
        d = self.num_states
        
        try: 
            #Assert correct shape and datatype
            if ham is not None:
                assert(ham.getHamMatrix().shape == (d,d))
            else:
                assert(self.H.getHamMatrix().shape == (d,d))
        except:
            print("Required shape: ", (d,d))
            raise exception_hamiltonian_compatibility

    """GET methods"""

    def getInitializedHamiltonian(self) -> bool:
        return self.initialized_H
    
    def getInitializedLattice(self) -> bool:
        return self.initialized_lat
    
    def getParameters(self) -> ParameterSet:
        return self.parameters
    
    #def getG(self) -> greenstensor_types:
     #   return self.G
    
    def getLattice(self) -> Lattice:
        self.assertInitializedLattice()
        return self.lat
    
    def getHamiltonian(self) -> Hamiltonian:
        self.assertInitializedHamiltonian()
        return self.H

    #Moving herm metric to StateCollection class
    #def getHermiticityMetric(self) -> ndarray:
    #    return self.P
    
    def getEigenStates(self) -> StateCollection:
        self.assertInitializedHamiltonian()
        assert(self.eigen_states is not None)
        
        return self.eigen_states
    
    def getNumStates(self) -> dtypes_int:
        return self.num_states

    """SET methods"""

    #def setInitializedHamiltonian(self, ini : bool) -> None:
    #    self.initialized_H = ini
    #These two bools should be private, right?
    #def setInitializedLattice(self, ini : bool) -> None:
    #    self.initialized_lat = ini

    def setParameters(self, parameters : ParameterSet) -> None:
        self.assertParametersStandard()
        self.parameters = parameters

    #def setG

    def setLattice(self, lattice : Lattice) -> None:
        self.lat = lattice
        self.initialized_lat = True
    
    def setHamiltonian(self, hamiltonian : Hamiltonian) -> None:
        self.assertHamMatrixProperties(hamiltonian)
        self.H = hamiltonian
        self.initialized_H = True

    #Moving herm metric to StateCollection class
    #def setHermiticityMetric(self, P : ndarray) -> None:
    #    self.assertInitializedHamiltonian()
    #    self.P = P

    def setEigenStates(self, states : StateCollection) -> None:
        self.assertInitializedHamiltonian()
        self.eigen_states = states

    def setNumStates(self, parameters : ParameterSet = None, N : dtypes_int = None, num_excitations : dtypes_int = None) -> None:
        """
        Set number of states of Atom Array. Either supply ParameterSet or N and num_excitations.

        --- INPUT 

         - parameters (ParameterSet)
         - N (dtypes_int) - see utils for Literal of ints
         - num_excitations (dtypes_int)
        """
        if parameters and N and num_excitations is None:
            self.num_states = number_states(self.parameters["N"][0], self.parameters["num_excitations"][0])
        elif N and num_excitations is not None:
            self.num_states = number_states(N, num_excitations)
        elif parameters is not None:
            self.num_states = number_states(parameters["N"][0], parameters["num_excitations"][0])
        else:
            print("When trying to set number of states: ")
            raise exception_required_parameters

    """Functions"""

    def initEigenStates(self, sort : bool = False) -> None:
        """
        Initialize eigenstates. To be called after initializing Hamiltonian. 

        Saves the right eigenvectors, default output of Hamiltonian eigendecomp, in the states dictionary. 
        """

        self.assertInitializedHamiltonian()

        N = self.num_states

        if not self.H.isDecomposed():
            self.H.eigenDecomposition()

        val, vec = self.H.getEigenDecomp()
        if sort:
            self.H.getDecayRates(sort=sort)
            sorted_index = self.H.getSortedIndex()
        else:
            sorted_index = range(N)

        #init hermiticity metric: #moving herm metric to StateCollection class. 
        #self.P = hermiticity_metric(vec)

        #Extract right eigenvectors from decomposition, which are the columns of vec.
        eigen_states = StateCollection({sorted_index[j]: State(N, sorted_index[j], val[sorted_index[j]], vec[:,sorted_index[j]]) for j in range(N)})

        #Manually set hermiticity metric in StateCollection, because init method is introducing too many errors, see workbench.ipynb
        eigen_states.setHermiticityMetric(hermiticity_metric(vec))

        self.eigen_states = eigen_states

    def updateEigenStates(self, new_states : StateCollection, **kwargs) -> None:
        """
        Update the eigen states during e.g. parameter sweep

        See StateCollection.updateStates() for a list of kwargs
        """
        
        self.eigen_states.updateStates(new_states.getStates(), **kwargs)