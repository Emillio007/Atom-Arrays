"""
@Author: tzs820, s240670, Emil Henningsen.

Utilities used throughout the subradiant states program.
"""

#Standard lib imports:
import warnings
from copy import deepcopy
from warnings import warn
from typing import Literal, Union, TypedDict, get_args, Self, Type, List
from collections.abc import Callable
from functools import singledispatchmethod
from textwrap import wrap

#from __future__ import annotations

#Third party lib imports:
import numpy as np
from numpy import (
    int64, float64, complex128, 
    dtype, ndarray, 
    zeros, concatenate, linspace, array, argmax, argmin, where, 
    real, imag, 
    exp)
from numpy.linalg import inv, outer
from scipy.constants import pi
from scipy.special import comb


"""Static elements: """
ex = np.array([1, 0, 0])                                         #x-axis unit vector
ey = np.array([0, 1, 0])                                         #y-axis unit vector
ez = np.array([0, 0, 1])                                         #z-axis unit vector

"""Program elements: """
#Datatypes currently allowed for values:
dtypes_int = Union[int, int64]
dtypes = Union[int, float, int64, float64, complex128]

#Float precision tolerance:
tolerance = 1e-1

#Allowed lattices, hamiltonians and greens tensors:
lattice_types = Literal["linear", "broken", "circular"]
hamiltonian_types = Literal["block", "scalar", "full"]
greenstensor_types = Literal["free-space"]

#Parameter settings
paramspace_types = Literal["linear", "sigmoid", "single"]

#Lattice-specific settings:
circle_measures = Literal["inter", "radius"]
circle_polarizations = Literal["inwards", "radial alternating", "outwards", "azimuthal", "azimuthal alternating", "other"]

#State reorder methods:
reorder_methods = Literal["overlap", "traveling", "shortest"]

"""Warnings: """
warn_herm_metric_empty = "No hermiticity metric passed to function or ndarray is empty"
warning_pola_none = Warning("No polarization is provided, proceeds with ez for all dipoles. ")
warning_hamiltonian_eigendecomp = Warning("Hamiltonian has not been eigendecomposed, proceeding to do so.")
warning_hamiltonian_probabilities = Warning("Probabilities exceed 1.")
warning_supplied_dtype = Warning("The supplied p_start and p_end do not match data types or is wrong data type. Note: the allowed data types are np.float64 and np.complex128")
warning_space_type = Warning("Type of parameter space incompatible. Defaults to linear.")
warning_timed_paramspace_excess = Warning("Number of param steps exceeds the final value of time steps. The excess param steps are excluded.")

exception_not_yet_implemented = Exception("The requested feature has not been implemented yet. ")
e = Exception("Hamiltonian is empty. ")
exception_standard_parameters = Exception("Standard parameters, N and num_excitations, are not specified. ")
exception_required_parameters = Exception("The required parameters for the chosen type of Atom Array are not specified. ")
exception_hamiltonian_compatibility = Exception("The Hamiltonian is of incorrect shape and / or datatype. ")
exception_datatype = Exception("The provided datatype is not supported. ")
exception_timed_paramspace = Exception("Failed to create timed parameter space")

"""Functions: """
def cos(angle):
    return np.cos(angle)
def sin(angle):
    return np.sin(angle)

def sigmoid(x, x0, a):
    return 1/(1 + np.exp(-a*(x-x0)))

def dagger(arr : ndarray) -> ndarray:
    return np.conjugate(np.transpose(arr))

def inner_nh(one : ndarray, two : ndarray, P : ndarray = None, right_eigenvectors : bool = True) -> ndarray:
    """
    Compute the Hilbert-space inner product of two state vectors

    --- INPUT

    one - NDarray of shape (N,) or (N,N)
    two - NDarray of shape (N,) or (N,N)
    right_eigenvectors - bool (optional). True if both vector arrays passed are right eigenvectors. 
    
    --- OUTPUT

    out - Complex-valued scalar or NDarray of shape (N,N)
    """
    if right_eigenvectors and P is None:
        warnings.warn(warn_herm_metric_empty)
    elif right_eigenvectors and P is not None:
        #The NH inner product defined from hermiticity metric, P
        return dagger(one) @ P @ two
    else:
        #The usual hermitian inner product
        return dagger(one) @ two

def hermiticity_metric(vecs : ndarray) -> ndarray:
    """
    Construct hermiticity metric, P, from right eigenvectors (Check HFT article for definition)
    """

    vecsDagger = dagger(vecs)

    """Construct hermiticity correction matrix 'G' (we will call it P)"""
    P = inv(vecs @ vecsDagger)
    
    return P

def number_states(N : int, num_excitations : int) -> int:
    return comb(N, num_excitations, exact=True)