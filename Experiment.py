"""
Experiment class to run atom array experiments with dynamic and static parameters.

Auhtor: Emil Henningsen, tzs820, s240670
"""

from State import *
from AtomArray import *
from ParameterSet import *
from ExperimentData import *

from utils import *

class Experiment:

    """variables:"""

    atoms = None
    parameters = None

    num_states = None

    evolve_decay = None

    T_run = None
    T_steps = None
    T_steps_it = None

    #T_begin = []

    dtau = None
    taus = None

    """save data:"""

    integrated_eigenvalues = None
    data = None

    def __init__(self, atoms : AtomArray, parameters : ParameterSet, evolve_decay : bool, 
                 T_run : int, T_steps : int, time_parameters : bool = True):

        self.atoms = atoms
        self.parameters = parameters

        self.num_states = number_states(parameters["N"][0], parameters["num_excitations"][0])
        self.evolve_decay = evolve_decay
        self.T_run = T_run 
        self.T_steps = T_steps
        self.T_steps_it = range(T_steps)
        #self.T_begin = T_begin

        self.dtau = T_run / T_steps
        self.taus = linspace(0, T_run, T_steps)

        """ what data to save? """

        #The complete data dictionary:
        self.data = ExperimentData(self.T_steps, self.taus, parameters) 

        #init params:
        if time_parameters:
            self.timeParameters()
    
    def timeParameters(self):
        """ initialization instruction """

        #cast parameterspaces into the timed version, so they fit into the iteration cycle
        #set up iterators
        for param in self.parameters:
            parameter = self.parameters[param]
            if type(parameter) == Parameter:
                try:
                    parameter.createTimedSpace(self.T_steps)
                    iter(parameter)
                except:
                    print("Tried to create timed parameter space for parameter: ", parameter.getName())
                    raise

    def resetParameters(self):
        for param in self.parameters:
            if type(self.parameters[param]) == Parameter:
                iter(self.parameters[param])

    def doExperiment(self, trial_state_collection : StateCollection, sort_new_states_lowest : bool = False,
                     evolution_basis : Literal["natural", "adiabatic"] = "natural", euler_type : Literal["end", "mid", "RK4"] = "end",
                     ordering_method : Literal["traveling_then_shortest", "traveling", "shortest"] = "traveling_then_shortest",
                     eigen_compute_left : bool = True, eigen_compute_hermiticity : bool = True, 
                     trial_compute_left : bool = False, return_M : bool = False, return_diffs : bool = False) -> ExperimentData:
        
        #Reset parameters:
        self.resetParameters()

        atoms = self.atoms

        if euler_type == "RK4":
            return_M = True
            return_diffs = True

        #TODO: acquire lattice and hamiltonian from atomarray and implement different tracking methods

        self.integrated_eigenvalues = {key: 0j for key in list(atoms.getEigenStates().getStates().keys())}

        for time_stamp in self.T_steps_it:

            """ --- one complete iteration of update cycle --- """

            #if feeding lattice and hamiltonian specifically, update them here:
            
            #update AtomArray:
            atoms.updateLattice(lattice="linear")
            atoms.updateHamiltonian(hamiltonian="block")

            H = atoms.getHamiltonian()

            """ diagonalize new hamiltonian """

            H.eigenDecomposition()

            eigvals, eigvecs = H.getEigenDecomp() #updated eigenvectors for 
            
            """ StateCollection to update position of eigenstates """
            #eigenstates don't change for now, since we don't update hamiltonian in any way.
            #but let's still work it through and see if something goes wrong. 
            if sort_new_states_lowest:
                H.getDecayRates()
                indicies_sorted = H.getSortedIndex()
            else:
                indicies_sorted = range(self.num_states)

            new_states = StateCollection({j: State(self.num_states, j, eigvals[j], eigvecs[:,j]) for j in indicies_sorted})

            #good for slow evolution (adiabatic)
            match(ordering_method):
                case "traveling_then_shortest":
                    if time_stamp == 0:
                        atoms.updateEigenStates(new_states, method="traveling", 
                                                compute_hermiticity_metric = eigen_compute_hermiticity, compute_left_set = eigen_compute_left)
                    else:
                        atoms.updateEigenStates(new_states, method="shortest", 
                                                compute_hermiticity_metric = eigen_compute_hermiticity, compute_left_set = eigen_compute_left)
                case "traveling":
                    atoms.updateEigenStates(new_states, method="traveling", 
                                            compute_hermiticity_metric = eigen_compute_hermiticity, compute_left_set = eigen_compute_left)
                case "shortest":
                    atoms.updateEigenStates(new_states, method="shortest", 
                                            compute_hermiticity_metric = eigen_compute_hermiticity, compute_left_set = eigen_compute_left)

            #StateCollection with reordered Eigenstates.
            sc = atoms.getEigenStates()

            """check for -1 fase factors on eigenvectors and correct if any"""

            if time_stamp != 0:
                prev_right_vectors = self.data.getRightVectors(time_stamp-1)
                for i in range(self.num_states):
                    sign = np.sign(dagger(prev_right_vectors[i]) @ sc[i].getRightVector()).real
                    #print("at time ", time_stamp, " the sign of state ", i, " is ", sign)
                    if time_stamp != 1:
                        sc[i].signFactor(sign)
                    else:
                        sc[i].signFactor(sign, evolved=False)

            #ensure also the left vectors and hermiticity metric is computed after any corrections
            sc.initHermiticityMetric()
            sc.computeLeftSet()

            """ update loop 1 , evolve eigenstates """

            for i in range(self.num_states):

                """ add to total eigenvalues (i.e. integration) """

                #shift = sc[i].getEnergyShift()
                #decay = sc[i].getDecayRate()
                #val = shift - 1j*decay                                         # since lambda_n = J_n - Gamma_n / 2 , we need val_n = Re - 2 * Im ... NO WE DON'T
                val = sc[i].getValue()
                self.integrated_eigenvalues[i] += val * self.dtau

                """ evolve eigenstates """

                sc[i].evolveState(self.integrated_eigenvalues[i], evolve_decay=self.evolve_decay)

            """ evolve other states (Euler iteration of Schrodingers equation) in basis of eigenstates.
            For adiabatic basis evolution: by using time-deriva of right vectors, we get the updated metric for the left vectors.
            See notes p. 43 and following meeting 15/1"""

            match(evolution_basis):

                case "natural": #skal være NATURLIGE basis

                    match(euler_type):

                        case "end":
                            
                            for trial_state in trial_state_collection:

                                rvec = trial_state.getEvolvedRightVector()

                                rvec += -1j * H.getHamMatrix() @ rvec * self.dtau

                                trial_state.setEvolvedRightVector(rvec, compute_left=trial_compute_left)

                        #case "mid":
                            
                            """
                            if time_stamp != 0:
                                
                                for trial_state in trial_state_collection:
                                    
                                    prev_rvec = self.data[time_stamp-1]["trialstate_collection"][trial_state.getName()].getEvolvedRightVector()
                                    rvec = trial_state.getEvolvedRightVector()

                                    rvec += -1j * """

                            

                        case _:
                            print("when trying Euler iteration type: ", euler_type, " with evo basis type: ", evolution_basis)
                            raise exception_not_yet_implemented

                case "adiabatic":

                    match(euler_type):
                        
                        case "end": #see case "mid" for comments on operations
                            
                            if time_stamp != 0:
                                #previous_right_vectors = self.data.getEvolvedRightVectors(time_stamp-1) #list
                                previous_right_vectors = self.data.getRightVectors(time_stamp-1, arr = True) #trying without evolution for coupling product!
                            else:
                                previous_right_vectors = [sc[i].getRightVector() for i in range(self.num_states)]

#                            right_vectors_difference = [(previous_right_vectors[i] - sc[i].getEvolvedRightVector())/self.dtau for i in range(self.num_states)]
                            right_vectors_difference = [(previous_right_vectors[i] - sc[i].getRightVector())/self.dtau for i in range(self.num_states)]

                            if time_stamp != 0:

                                #previous_eigenstates = self.data.getEigenStateCollection(time_stamp-1)
                                #previous_integrated_eigenvalues = self.data[time_stamp-1]["integrated_eigenvalues"]
                                
                                for trial_state in trial_state_collection:

                                    #previous_trial_state = self.data[time_stamp-1]["trialstate_collection"][trial_state.getName()].getEvolvedRightVector()
                                    
                                    new_rvec = []
                                    for coefficient, m in zip(trial_state.getEvolvedRightVector(), range(self.num_states)): 

                                        #Prøver med ikke-integrerede egenværdier ... Vi udfører jo nemlig en integration her, så vi skal vel ikke dobbelt?
                                        eigvel = sc[m].getValue()

                                        #End-point euler it (not integrated eigvals)
                                        first_term = -1j * (eigvel*coefficient) * self.dtau

                                        second_term = 0
                                        for n in range(self.num_states):

                                            #End-point euler it:
                                            #second_term += trial_state.getEvolvedRightVector()[n] * (sc[m].getEvolvedLeftVector() @ right_vectors_difference[n]) * self.dtau
                                            second_term += trial_state.getEvolvedRightVector()[n] * (sc[m].getLeftVector() @ right_vectors_difference[n]) * self.dtau

                                        coefficient += first_term - second_term
                                        new_rvec.append(coefficient)

                                    trial_state.setEvolvedRightVector(array(new_rvec), compute_left=trial_compute_left)

                        case "mid":

                            #for first time step, previous iteration vectors are the same, therefore difference term drops out
                            if time_stamp not in [0, 1]:
                                
                                previous_right_vectors = self.data.getRightVectors(time_stamp-1) #list
                                previous_right_vectors2 = self.data.getRightVectors(time_stamp-2)
                            
                            else:
                                
                                previous_right_vectors = [sc[i].getRightVector() for i in range(self.num_states)]
                                previous_right_vectors2 = [sc[i].getRightVector() for i in range(self.num_states)]

                            #first order differentiation of kets
                            right_vectors_difference = [(previous_right_vectors[i] - sc[i].getRightVector())/self.dtau for i in range(self.num_states)] 
                            right_vectors_difference2 = [(previous_right_vectors2[i] - previous_right_vectors[i])/self.dtau for i in range(self.num_states)]

                            if return_M:
                                #Test matrix to check if negative imaginary eigenvalues:
                                M = zeros((self.num_states, self.num_states), dtype=complex128)

                            #now, let's try mid point Euler iteration instead of endpoint, i.e. skip first two iterations to have a mid point (with diff vectors)
                            if time_stamp != 0 and 1:
                                
                                previous_eigenstates = self.data.getEigenStateCollection(time_stamp-1)
                                
                                #update coefficients ... this can probably be done faster in numpy arrays 
                                for trial_state in trial_state_collection:
                                    
                                    previous_trial_state = self.data[time_stamp-1]["trialstate_collection"][trial_state.getName()].getEvolvedRightVector()
                                    
                                    new_rvec = []
                                    for coefficient, m in zip(trial_state.getEvolvedRightVector(), range(self.num_states)): 
                                        
                                        prev_coefficient = previous_trial_state[m]

                                        #Prøver med ikke-integrerede egenværdier ... Vi udfører jo nemlig en integration her, så vi skal vel ikke dobbelt?
                                        eigvel = sc[m].getValue()
                                        prev_eigval = previous_eigenstates[m].getValue()
                                        
                                        #mid-point euler it (not integrated eigvals)
                                        first_term = -1j * (prev_eigval*prev_coefficient + eigvel*coefficient)/2 * self.dtau

                                        second_term = 0
                                        for n in range(self.num_states):
                                            
                                            #mid-point euler it with correct weighted average
                                            second_term += (previous_trial_state[n] * previous_eigenstates[m].getEvolvedLeftVector() @ right_vectors_difference2[n] + \
                                                            trial_state.getEvolvedRightVector()[n] * sc[m].getEvolvedLeftVector() @ right_vectors_difference[n]) / 2 * self.dtau
                                            if return_M:
                                                if m != n:
                                                    M[m,n] -= second_term
                                        
                                        coefficient += first_term - second_term
                                        
                                        if return_M:
                                            M[m,m] = coefficient
    
                                        new_rvec.append(coefficient)
    
                                    trial_state.setEvolvedRightVector(array(new_rvec), compute_left=trial_compute_left)


                            #Check M
                            #decay_rates_M = -2*imag(np.linalg.eigvals(M))
                            #print("at time: ", time_stamp, " the negative decay rates are: ", np.where(decay_rates_M<0))
                        
                        case "RK4":
                            
                            if time_stamp not in [0, 1]:
                                previous_right_vectors = self.data.getRightVectors(time_stamp-1, arr=True)
                                previous_right_vectors_2 = self.data.getRightVectors(time_stamp-2, arr=True)
                            else:
                                previous_right_vectors = array([sc[i].getRightVector() for i in range(self.num_states)])
                                previous_right_vectors_2 = array([sc[i].getRightVector() for i in range(self.num_states)])


                            right_vectors_difference = ( previous_right_vectors - array([sc[i].getRightVector() for i in range(self.num_states)]) ) / self.dtau
                            

                            if time_stamp not in [0, 1]:

                                previous_eigenstates = self.data.getEigenStateCollection(time_stamp-1)
                                #get right_vectors_difference_2 from data

                                for trial_state in trial_state_collection:

                                    previous_trial_state = self.data[time_stamp-1]["trialstate_collection"][trial_state.getName()].getEvolvedRightVector()
                                    
                                    # do RK4 iteration of trial state

                                    trial_state.setEvolvedRightVector(array(new_rvec), compute_left=trial_compute_left)


                        case _:
                            print("when trying Euler iteration type: ", euler_type, " with evo basis type: ", evolution_basis)
                            raise exception_not_yet_implemented

                case _:

                    print("when trying evolution basis: ", evolution_basis)
                    raise exception_not_yet_implemented

            """ save self.data """

            #Remember to use deepcopies, as otherwise the new self.data is written in same instance overwriting the old self.data. 
            self.data[time_stamp]["integrated_eigenvalues"] = deepcopy(self.integrated_eigenvalues)
            self.data[time_stamp]["eigenstate_collection"] = deepcopy(sc)
            self.data[time_stamp]["trialstate_collection"] = deepcopy(trial_state_collection)
            if return_M:
                self.data[time_stamp]["full_matrix"] = deepcopy(M)
            if return_diffs:
                self.data[time_stamp]["differentials"] = deepcopy(right_vectors_difference)

            #For large values of N and T_steps, the algorithm should regularly write to file, as e.g. N=50 and steps=5000 requires 6,5 GB of memory otherwise.
            #This is one solution. Otherwise, try to set up program on HPC. 

            """ Update parameters for next iteration """

            for param in self.parameters:
                parameter = self.parameters[param]
                if type(parameter) == Parameter:
                    next(parameter)

        return self.data