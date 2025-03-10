"""
The ParameterSet class for use in Experiments with Atom Array.

Author: Emil Henningsen, tzs820, s240670
"""

from utils import *
from Parameter import *

class ParameterSet(TypedDict, total=False):

    #Basic system parameters:
    N: Parameter
    num_excitations: Parameter
    greenstensor_type: greenstensor_types
    
    #Parameters for linear lattice (unidirectional polarizations)
    distance: Parameter
    direction: Parameter
    angle: Parameter
    polarization: Parameter #Maybe only one of the two??

    #For broken linear
    broken_angle : Parameter

    #Parameters for circular lattice
    circle_measure: circle_measures
    circle_polarization: circle_polarizations
    radius: Parameter


