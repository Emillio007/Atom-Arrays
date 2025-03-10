"""
The Parameter class for use in Experiments with Atom Array.

Author: Emil Henningsen, tzs820, s240670

TODO:

 - Fix datatype problem (must handle both np.array dtypes and type(x) types.)... half-way done
"""

from utils import *

class Parameter:

    """Internal storage:"""
    name = None
    static = False

    p0, p1 = None, None
    datatype = None

    paramspace_type = None
    paramspace = None           #Must be saved as numpy array
    steps = 0
    steps_it = None

    # for use when iterating over a time-interval:
    timed = False
    T_begin = None
    timed_paramspace = None
    timed_steps = None

    """private variables"""
    
    a = None

    """ init """

    def __init__(self, 
                 name : str, 
                 static : bool,
                 p_start : dtypes, 
                 p_end : dtypes = None, 
                 steps : int = None, 
                 space : paramspace_types = "single", 
                 T_begin : int = 0):
        """Module for using parameters in AtomArray. So far, only 1D parameter spaces are supported.
        
        --- INPUT

         - name (str) - Name of this parameter, e.g. angle
         - static (bool) - Toggle for static parameter, in this case, no p_end or steps is to be supplied
         - p_start (float64 or complex128) - starting point, e.g. 0
         - p_end (float64 or complex128) - end point, e.g. pi/2
         - steps (int) - number of steps
         - space (str, literal) - Type of parameter space: ["linear", "sigmoid", "single"]. Default: "single"

        --- OUTPUT

        Initialized Parameter object.
        """

        self.name = name
        self.static = static
        self.paramspace_type = space
        self.T_begin = T_begin

        if static:
            self.assertDataType(p_start)
            self.p0 = p_start
            self.datatype = type(p_start)
            self.createSingle()
        else:
            self.steps = steps
            self.steps_it = range(steps)
            #Assert correct datatypes and store
            self.assertDataType(p_start)
            self.assertDataType(p_end)

            if type(p_start) == type(p_end):
                self.p0 = p_start
                self.p1 = p_end
                self.datatype = type(p_start)
            else:
                warn(warning_supplied_dtype)

            #Fill out parameterspace:
            match space:
                case "linear":
                    self.createLinearSpace()
                case "sigmoid":
                    self.createSigmoidSpace()
                case _:
                    self.createLinearSpace()
                    self.paramspace_type = "linear"

    """Overriding operators: """

    def __getitem__(self, key):
        """ if timed, returns the timed_paramspace[key] """
        if key == "auto":
            key = self.getCurrentIteration()
        if self.timed: 
            return self.timed_paramspace[key]
        else:
            return self.paramspace[key]
    
    def __setitem__(self, key, value):
        """ if timed, returns the timed_paramspace[key] """
        if key == "auto":
            key = self.getCurrentIteration()
        if self.timed: 
            self.timed_paramspace[key] = value
        else:
            self.paramspace[key] = value

    def __iter__(self):
        """ if timed, iterates the timed_paramspace[key] """
        self.a = 0
        return self
    
    def __next__(self):
        """ if timed, iterates the timed_paramspace[key] """
        if self.timed:
            if self.a < self.timed_steps:
                val = self.timed_paramspace[self.a]
                self.a += 1
                return val
            else:
                raise StopIteration
        else:
            if self.a < self.steps:
                val = self.paramspace[self.a]
                self.a += 1
                return val
            else:
                raise StopIteration

    """Assertions"""

    def assertDataType(self, x : dtypes) -> None:
        """
        Assert if x is of allowed datatypes. See utils for allowed datatypes.
        """        
        if type(x) in get_args(dtypes):
            pass
        elif type(x) == ndarray:
            assert(x.dtype in get_args(dtypes))
        else:
            print(type(x), x)
            raise exception_datatype

    """ helper functions """

    def createSingle(self) -> None:
        self.paramspace = array([self.p0])

    def createLinearSpace(self) -> None:
        """
        Create linear space (from Numpy's linspace function) and store internally
        """
        self.paramspace = linspace(self.p0, self.p1, self.steps)#, dtype=self.datatype)

    def createSigmoidSpace(self) -> None:
        """
        Create sigmoid space (See utility function sigmoid) and store internally.
        Sigmoid from utils runs from 0 to 1, so p0 is added for zero-point and p1 is factored for scaling, 
        i.e. now sigmoid runs from p0 to p1.
        The sigmoid x0 value is halfway through the steps, i.e. it is centered around n/2.
        The sigmoid a value is 4*pi/n, which iteratively was shown to accelerate smoothly.
        """
        x0 = self.p0    #Starting parameter value
        x1 = self.p1    #Ending parameter value
        n = self.steps          #Number of steps
        n_it = self.steps_it    #Iterator of steps
        if x1 < x0: #if starting from greater value to lower
            self.paramspace = x1 + sigmoid(array(list(reversed(n_it))), n/2, (4*pi)/n) * (x0-x1)
        else:   #regular direction:
            self.paramspace = x0 + sigmoid(array(n_it), n/2, (4*pi)/n) * (x1-x0) #Start at x0, i.e. 0->x0, and rescale sigmoid 0-1 to 0-x1

    def createTimedSpace(self, T_steps : int, T_begin : int = None) -> None:
        """
        Insert parameter space in a timed setting, i.e. when using in experiment. 
        Stores a timed_parameterspace, which can be used to 

        --- INPUT

         - T_steps (int) - the number of time steps in total
         - T_begin (int) - When to begin the parameterspace changing - defaul=0
        
        """

        if T_begin is not None:
            self.T_begin = T_begin

        timed_space = self.paramspace

        begin_param_value = self.paramspace[0]
        end_param_value = self.paramspace[-1]

        # if T_begin is 0, no need for a front piece

        if self.T_begin == 0:
            front_piece = None
        else:
            front_piece = array([begin_param_value for time in range(self.T_begin)])
        
        # if number of param steps exceeds number of time steps, don't make backpiece and cut off excess

        if self.T_begin+self.steps >= T_steps:
            back_piece = None
            timed_space = timed_space[:T_steps-self.T_begin]
            warnings.warn(warning_timed_paramspace_excess)
        else:
            back_piece = array([end_param_value for time in range(self.T_begin+self.steps, T_steps)])

        # check if front and/or back is provided and concatenate, otherwise return 

        if front_piece is None and back_piece is None:
            pass
        elif front_piece is None:
            timed_space = concatenate((timed_space, back_piece))
        elif back_piece is None:
            timed_space = concatenate((front_piece, timed_space))
        else:
            timed_space = concatenate((front_piece, timed_space, back_piece))

        #Store the new timed paramspace
        self.timed_steps = T_steps
        self.timed_paramspace = timed_space

        #finally raise the timed flag:
        self.timed = True

    """GET methods"""

    def getName(self) -> str:
        return self.name
    
    def getStart(self) -> dtypes:
        return self.p0
    
    def getEnd(self) -> dtypes:
        return self.p1
    
    def getSteps(self) -> int:
        return self.steps
    
    def getDataType(self) -> dtype:
        return self.datatype
    
    def getSpaceType(self) -> str:
        return self.paramspace_type
    
    def getSpace(self) -> ndarray:
        return self.paramspace
    
    def getStatic(self) -> bool:
        return self.static
    
    def getCurrentIteration(self) -> int:
        if self.a is None:
            return 0
        else:
            return self.a

    """SET methods"""

    def setName(self, name : str) -> None:
        self.name = name
    
    def setStart(self, p_start : dtypes) -> None:
        self.p0 = p_start
    
    def setEnd(self, p_end : dtypes) -> None:
        self.p1 = p_end
    
    def setSteps(self, steps : int) -> None:
        self.steps = steps
        self.steps_it = range(steps)
    
    def setDataType(self, dtype : dtypes) -> None:
        self.datatype = dtype
    
    def setSpaceType(self, paramspace_type : paramspace_types) -> None:
        self.paramspace_type = paramspace_type
    
    def setSpace(self, space : ndarray) -> None:
        self.paramspace = space
    
    def setStatic(self, static : bool) -> None:
        self.static = static
 
    def setCurrentIteration(self, i : int) -> None:
        self.a = i