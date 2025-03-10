"""
The ExperimentData class - controlled dict structure to store data from Experiment

Author: Emil Henningsen, tzs820, s240670

"""

from utils import *
from StateCollection import StateCollection
from ParameterSet import ParameterSet

class ExperimentData:

    _T_steps = None

    _data = {}
    _parameters = None

    def __init__(self, T_steps : int, taus : ndarray, parameters : ParameterSet):
        self._T_steps = T_steps
        self._data = {
            time_stamp: {
                "time": taus[time_stamp],               #Preentered time corresponding to current time_stamp
                "integrated_eigenvalues": dict,         #Save integrated eigenvalues of current time_stamp here
                "differentials": ndarray,                      #Save couplings of current time_stamp here
                "full_matrix": ndarray,                         #Save full evolution matrix of time_stamp here
                "eigenstate_collection": StateCollection,     #Save EigenStateCollection of current time_stamp here
                "trialstate_collection": StateCollection
            } for time_stamp in range(T_steps)
        }
        self._parameters = parameters

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value

    def getNum(self) -> int:
        #For now just default to N as number of states
        return self._parameters["N"][0]

    """ GET """

    def getTimes(self, arr : bool = True) -> Union[list[float64], ndarray]:
        t = [self._data[time]["time"] for time in range(self._T_steps)]
        if arr:
            t = array(t)
        return t

    """ variations for eigenvalues """

    def getEigenvalue(self, time_stamp : int, state : int) -> dtypes:
        return self._data[time_stamp]["eigenstate_collection"][state].getValue()

    def getEigenvalues(self, time_stamp : int, arr : bool = True) -> Union[list, ndarray]:
        """ if arr, returns ndarray otherwise list of complex eigenvalues at time stamp """
        N = self.getNum()
        vals = [self.getEigenvalue(time_stamp, state) for state in range(N)]
        if arr:
            vals = array(vals)
        return vals

    def getIntegratedEigenvalues(self, time_stamp : int) -> dict:
        """
        get eigenvalues sorted in dict such that 

        {key: complex_value(key) for key in list(atoms.getEigenStates().getStates().keys())}
        """
        return self._data[time_stamp]["integrated_eigenvalues"]
    
    #def getIntegratedEigenvaluesTotal(self) -> dict:
    #    return {time_stamp: }

    def getEigenvalueTotal(self, state : int, arr : bool = True) -> Union[list, ndarray]:
        """Return Eigenvalue for all times. If arr is True, is returned in numpy array with shape (T_steps)"""
        value = [self._data[time_stamp]["eigenstate_collection"][state].getValue() for time_stamp in range(self._T_steps)]
        if arr:
            value = array(value)
        return value

    def getEigenvaluesTotal(self, arr : bool = True) -> Union[list, ndarray]:
        """
        Return Eigenvalues for all times. if arr is True, is returned in numpy array with shape (T_steps, N)

        Otherwise is returned as list
        """
        N = self.getNum()
        values = [[self._data[time_stamp]["eigenstate_collection"][state].getValue() for state in range(N)] for time_stamp in range(self._T_steps)]
        if arr:
            values = array(values)
        return values

    """ variations for names """

    def getName(self, time_stamp : int, state : int) -> int:
        return self._data[time_stamp]["eigenstate_collection"][state].getName()

    def getNames(self, time_stamp : int) -> list:
        N = self.getNum()
        names = [self._data[time_stamp]["eigenstate_collection"][state].getName() for state in range(N)]
        return names

    """ variations for states """

    def getEigenState(self, time_stamp : int, state : int):
        return self._data[time_stamp]["eigenstate_collection"][state]
    
    def getEigenStateTotal(self, state : int) -> list:
        """get state at all times as list"""
        states = [self._data[time_stamp]["eigenstate_collection"][state] for time_stamp in range(self._T_steps)]
        return states

    def getEigenStateCollection(self, time_stamp : int):
        return self._data[time_stamp]["eigenstate_collection"]

    def getTrialState(self, time_stamp : int, state : int):
        return self._data[time_stamp]["trialstate_collection"][state]
    
    def getTrialStateTotal(self, state : int) -> list:
        """get state at all times as list"""
        states = [self._data[time_stamp]["trialstate_collection"][state] for time_stamp in range(self._T_steps)]
        return states
    
    def getTrialStateEvolvedVectorTotal(self, state : int, arr : bool = False, in_cols : bool = False) -> Union[list, ndarray]:
        vec = [self._data[time_stamp]["trialstate_collection"][state].getEvolvedRightVector() for time_stamp in range(self._T_steps)]
        if arr:
            vec = array(vec)
        if in_cols:
            if arr:
                vec.transpose()
                
            else:
                print("When trying getTrialStateEvolvedVectorTotal with in_cols, ", in_cols)
                raise exception_not_yet_implemented
        return vec

    def getTrialStateCollection(self, time_stamp : int):
        return self._data[time_stamp]["trialstate_collection"]

    """ variations for vectors 
    
        TODO: @singledispatch s.t. only one method to call, but varying parameters, e.g. no state should return all states then
        ... actually no need for singledispatch. Implement one method, getVector(time, state, space="right") and then do match cases and implement all possible combinations.
    """

    def getRightVector(self, time_stamp : int, state : int) -> ndarray:
        """ get right vector at time stamp """
        return self._data[time_stamp]["eigenstate_collection"][state].getRightVector()
    
    def getRightVectorTotal(self, state : int, arr : bool = True) -> list[ndarray]:
        """ get right vector at time stamp """
        vec = [self._data[time_stamp]["eigenstate_collection"][state].getRightVector() for time_stamp in range(self._T_steps)]
        if arr:
            vec = array(vec)        
        return vec

    def getLeftVector(self, time_stamp : int, state : int) -> ndarray:
        """ get left vector at time stamp """
        return self._data[time_stamp]["eigenstate_collection"][state].getLeftVector()
    
    def getRightVectors(self, time_stamp : int, arr : bool = False) -> Union[list[ndarray], ndarray]:
        """ get right vector at time stamp """
        x = [self._data[time_stamp]["eigenstate_collection"][state].getRightVector() for state in range(self.getNum())]
        if arr:
            x = array(x)
        return x

    def getRightVectorsAllTimes(self, arr : bool = True) -> Union[list[ndarray], ndarray]:
        """Get all right vectors at all times (T_steps, num_states, N)"""
        x = [[self._data[time_stamp]["eigenstate_collection"][n].getRightVector() for n in range(self.getNum())] for time_stamp in range(self._T_steps)]
        if arr:
            x = array(x)
        return x

    def getLeftVectors(self, time_stamp : int) -> list[ndarray]:
        """ get left vector at time stamp """
        return [self._data[time_stamp]["eigenstate_collection"][state].getLeftVector() for state in range(self.getNum())]

    def getLeftVectorsAllTimes(self, arr : bool = True) -> Union[list[ndarray], ndarray]:
        """Get all left vectors at all times (T_steps, num_states, N)"""
        x = [[self._data[time_stamp]["eigenstate_collection"][n].getLeftVector() for n in range(self.getNum())] for time_stamp in range(self._T_steps)]
        if arr:
            x = array(x)
        return x

    def getEvolvedRightVectors(self, time_stamp : int) -> list[ndarray]:
        return [self._data[time_stamp]["eigenstate_collection"][state].getEvolvedRightVector() for state in range(self.getNum())]

    def getHermiticityMetric(self, time_stamp : int) -> ndarray:
        return self._data[time_stamp]["eigenstate_collection"].getHermiticityMetric()

    """ SET """

    def setIntegratedEigenvalues(self, time_stamp : int, eigvals : dict):
        """
        set eigenvalues sorted in dict such that 

        {key: complex_value(key) for key in list(atoms.getEigenStates().getStates().keys())}
        """
        self._data[time_stamp]["integrated_eigenvalues"] = deepcopy(eigvals)
