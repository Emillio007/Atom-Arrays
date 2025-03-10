"""
The StateCollection class. Used to keep a collection of states, e.g. eigenstates, and update positions in dictionary upon parameter updates. 

Author: Emil Henningsen, tzs820, s240670
"""

from utils import *
from State import State

class StateCollection:

    _states : dict = None

    #Hermiticity Metric moved from AtomArray here:
    _P = None

    """constructor"""
    def __init__(self, states : dict):
        """
        Initialize StateCollection object. 

        --- INPUT

         - states (dict) - Dictionary of states in collection, keys must be numbers, e.g. the State names. 
        
        """

        self._states = states
        self.initHermiticityMetric()

    """overloading operators"""

    def __getitem__(self, key : int) -> State:
        return self._states[key]
    
    def __setitem__(self, key : int, value : State) -> None:
        self._states[key] = value

    def __iter__(self):
        """
        Create iterator of StateCollection. Iterates through N states as [0,N], like range(N).
        """
        self.a = 0
        return self
    
    def __next__(self):
        if self.a < len(self._states):
            state = self._states[self.a]
            self.a += 1
            return state
        else:
            raise StopIteration

    """GET methods"""
    def getStates(self) -> dict:
        """
        Get the Dictionary of States. keys in [0, num_states] and values State(name=key).
        """
        return self._states

    def getState(self, key : int) -> State:
        return self._states[key]

    def getHermiticityMetric(self) -> ndarray:
        return self._P

    """SET methods"""
    def setStates(self, states : dict) -> None:
        self._states = states

    def setState(self, key : int, value : State) -> None:
        self._states[key] = value

    def setHermiticityMetric(self, P : ndarray) -> None:
        self._P = P

    """Functions: """

    def initHermiticityMetric(self, method : Literal["product", "summation"] = "product") -> None:
        n = len(self._states[0].getVector())
        #Seems the summation way introduces too much numerical error? Probably the inversion routine introducing too much error. Update: Summation is invalid, as inversion does not distribute across sum.
        W = ndarray((n,n), dtype=complex128)
        
        if method == "summation":
            for state in self:
                W += state.getHermiticityMetricSingle()
        elif method == "product":
            vecs = []
            for state in self:
                vecs.append(state.getRightVector())
                #W += 
            vecs = array(vecs).transpose()
            W = hermiticity_metric(vecs)
        
        self._P = W

    def computeLeftSet(self, right_set : list = None, hermiticity_metric : ndarray = None):
        if right_set is None:
            right_set = self
        if hermiticity_metric is None:
            hermiticity_metric = self._P

        for state in right_set:
            state.setLeftVector(dagger(state.getRightVector()).transpose() @ hermiticity_metric)

    def updateStates(self, new_states : dict, index_list : list = None, method : reorder_methods = "traveling", 
                     fast : bool = False, compute_hermiticity_metric : bool = True, compute_left_set : bool = True) -> None:
        """
        Update states: reorder the states corresponding to new index list computed via chosen method and update to new states.
        The states will always stay at the same position in the dict of States, i.e. 0'th state is always at key=0 this way.

        --- INPUT
         - new_states (dict). Dictionary of States - new states to reorder.
         - index_list (list of ints) optional. The new indicies to arrange the states dictionary. If None, the list is computed using chosen method. 
         - method (Literal(str)) optional. The chosen method to reorder after. See utils for allowed methods. "overlap" is standard. 
         - fast (bool) optional. Toggle faster computation of overlap (compared with looping method). DefaulT: False. Not implemented yet!

        The "traveling" traveling salesman reorder method not yet implemented. This is (probably) particularly useful, when updating states whose eigenvalues lay nicely on a curve.
        """

        if index_list is None:
            match(method):
                case "traveling":
                    index_list = self.indexListByTraveling(new_states)
                case "shortest":
                    index_list = self.indexListByShortest(new_states)
                case "overlap":
                    index_list = self.indexListByOverlap(new_states, fast)
                case _:
                    print("When trying StateCollection reorder method: ", method)
                    raise exception_not_yet_implemented

        new_states = {j: new_states[i] for i, j in zip(index_list, range(len(index_list)))}

        self._states = new_states

        if compute_hermiticity_metric:
            self.initHermiticityMetric()
        if compute_hermiticity_metric and compute_left_set:
            self.computeLeftSet()
        elif compute_left_set and not compute_hermiticity_metric:
            print("When trying update states with compute_left_set, ", compute_left_set, " and compute_hermiticity_metric, ", compute_hermiticity_metric)
            raise exception_not_yet_implemented

    #Overlap computation:

    def indexListByOverlap(self, new_states : dict, fast : bool = False) -> list:
        """
        Compute list of indicies of new_states by overlap. I.e. new_index_i = max_j (bra(L_i)*ket(R_j))
        
        --- INPUT

         - new_states (dict). Dictionary of States - new states to compute overlap of.
         - fast (bool) optional. Toggle faster computation (to compare with looping method). DefaulT: False. Not implemented yet!
        
        Fast means doing matmul: ({bra(L_i)}) @ ({ket(R_j)}) instead of looping through each state and computing overlaps.
        """
    
        old_states = self._states
        new_indicies = []

        if fast:
            
            #See HFTpred_animated_v2 for overlap method using U @ V and then going through each column to find argmax
            #for now, not implemented

            raise exception_not_yet_implemented
        
        else:

            for i in old_states:
                
                overlap = []

                for j in new_states:
                    overlap.append(old_states[i] @ new_states[j])
                
                overlap = array(overlap)
                #print("overlaps of old_state ", i, " is ", overlap)

                new_indicies.append(int(argmax(overlap)))   #finding argmax of overlap and appending argument cast from np.int64 to int to use as key in new states.

        return new_indicies
    
    #Traveling salesman computation:

    def masker(self, i, n) -> list:
            return [True if v in i else False for v in range(n)]

    def distanceMatrix(self, values_a, values_b):
            n_points = len(values_a)
            assert(n_points == len(values_b))
            D = zeros((n_points, n_points))
            for i in range(n_points):
                for j in range(n_points):
                    D[i,j] = abs(values_a[i]-values_b[j])
            return D

    def indexListByTraveling(self, new_states : dict) -> list:
        """ TODO: description, for now, see workstation_2 """
    
        n_points = len(new_states)
        values = [new_states[i].getValue() for i in range(n_points)]

        #construct distance matrix
        D = self.distanceMatrix(values, values)

        stay_away_threshold = 10000.
        stay_away_matrix = np.full(n_points, stay_away_threshold)

        #search for best path:
        routes = {}

        for i in range(n_points):

            #route = np.zeros(n_points) #uncomment for NumPy only implementation.. but is slower
            route = []
            route_distance = np.float64(0)

            #We start at i'th index
            #route[0]= i    #uncomment for NumPy only implementation.. but is slower
            route.append(i)
            current_point = i

            for j in range(n_points-1): #we already have the first point
                
                next_points = D[:,current_point]
                
                #avoid going back to same points by masking them:
                mask = self.masker(route, n_points)
                
                #Go to next point, with the smallest distance:
                current_point = np.argmin(np.where(mask, stay_away_matrix, next_points))

                #route[j] = current_point   #uncomment for NumPy only implementation.. but is slower
                route.append(current_point)
                route_distance += next_points[current_point]
            
            #routes[i] = (route, route_distance)    #uncomment for NumPy only implementation.. but is slower
            routes[i] = (array(route), route_distance)

        route_distances = [routes[route][1] for route in routes]
        optimal_route = argmin(route_distances)
        ru = routes[optimal_route][0]

        return ru #return optimal route as new index list
    
    def indexListByShortest(self, new_states : dict):
        """ TODO: description, for now, see workstation_2 """
        old_states = self._states
        n_points = len(new_states)
        values = [new_states[i].getValue() for i in range(n_points)]
        values_before = [old_states[i].getValue() for i in range(n_points)]

        D = self.distanceMatrix(values, values_before)
        route = [] #just deliver one single route, because we are certain that it is the shortest
        for i in range(n_points):
            next_point = argmin(D[:,i])
            route.append(next_point)
        return array(route)