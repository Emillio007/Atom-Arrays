"""
Author: tzs820, s240670, Emil Henningsen.

GREEN'S TENSOR (in free space)
In this section, the Green's tensor for the i'th dipole w.r.t. the j'th dipole is calculated.
The assumptions include:
 * Markovian regime: I.e. the bandwidth of atomic transition, w0, is small and the retardation of the field is negligible.
 * Free space: I.e. the Green's tensor is calculated in free space, where the dipoles are tightly trapped, i.e. d < lambda_0.
 * The stroung coupling regime is avoided: I.e. the electromagnetic field is not strongly coupled to the atomic transition (CQED).
 * The dipoles are point-like particles: I.e. r_i and r_j can be treated as fixed position vectors. 
"""

from utils import *

def G_0(r, w0=1, dimensionless=True):
    """
    Green's tensor in free space. As given by eq. (6) of Asenjo-Garcia et al. 
    with option for dimensionless units. See bottom for definitions.

    ---Parameters:
    r: vector (3) of floats, displacement between the dipoles.
    w0 (optional): float, frequency of the atomic transition. Default: 1
    dimensionless (optional): bool, toggle dimensionless units. Default: True
    
    ---Return:
    G: array (3 x 3) of complex floats, Green's tensor in free space between two dipoles.

    ---

    Description:
    In the unitless case, no w needs to be provided. The variables are scaled such that:
        G = G' * k0
        r = r' * r0, r0 = 1/k0
    """
    from numpy.linalg import norm
    from numpy import e, outer, pi, identity
    from scipy.constants import c
    #from qutip import qeye
    
    k0 = w0/c
    r_norm = norm(r)
    if dimensionless:
        G = ((e**(1j * r_norm))/(4 * pi * r_norm**3)) * (
            (r_norm**2 + 1j * r_norm - 1) * identity(3) + 
            (-r_norm**2 -1j * 3 * r_norm + 3) * (outer(r,r)/r_norm**2)
            )
    else:
        G = ((e**(k0*r_norm*complex(0,1)))/(4 * pi * k0**2 * r_norm**3)) * (
            (k0**2 * r_norm**2 + k0 * r_norm * complex(0, 1) - 1)*identity(3) 
            + (-k0**2 * r_norm**2 - 3 * k0 * r_norm * complex(0, 1) + 3)*outer(r, r)/(r_norm**2)
            )
    return G

def fill_G(N, rij, w0=1, dimensionless=True, type : greenstensor_types = "free-space"):
    """
    Fill array of chosen Green's Tensor's for a given lattice with displacement vectors, rij. 

    ---Parameters:
    N: Integer number of dipoles.
    rij: array (N x N x 3) of floats, displacement vectors between dipoles in lattice.
    w0 (optional): float, dipole transition frequency. Default: 1
    dimensionless (optional): bool, default: True
    type (optional?): string, choose which Green's Tensor to compute. Default: "free-space"

    ---Return:
    G: array (N x N x 3 x 3) of complex floats, the Green's Tensor in free space for given lattice. 
    """
    from numpy import zeros, eye, pi

    G = zeros((N, N, 3, 3), dtype=complex)                        #Array of Green's tensors
    diagval = 1
    for i in range(N):                                            #Fill the array of Green's tensors
        for j in range(i, N):
            if i == j:
                G[i, j] = diagval*eye(3, dtype=complex)                   #Set diagonal values to 1 to avoid divergence.
            else:
                if type == "free-space":                          #Choose Green's Tensor by given optional argument.
                    G[i, j] = G_0(rij[i, j], w0, dimensionless)
                else:
                    print("Tried greens tensor type: ", type)
                    raise exception_not_yet_implemented

                G[j, i] = G[i, j]                                 #By construction, Green's tensor is symmetric.
    return G


def scalar(N, rij, w0=1, dimensionless=True):
    """
    G_zz scalar case of Green's Tensor function. 

    --- Parameters:

    N: integer, number of dipoles
    rij: array (N x N x 3) of floats, displacement vectors between dipoles.
    w0 (optional): float, transition frequency of dipoles. Default: True.
    dimensionless (optional): bool, compute dimensionless scalar case. Default: True.

    --- Return:
    
    h: array (N x N) of complex numbers, Green's function scalar. 

    ---

    Description:
    Below is the Scalar case, when using direction z as the only polarization direction of the dipoles. In that case, only the G_zz component is non-zero. 
    Furthermore, the function (denoted h) can be made dimensionless by introducing:
        r = r' * r0
        h = h' * k0
        r0 = 1/k0

    Reference: equation (6) of Asenjo-Garcia et al.

    """
    from numpy import zeros, e, pi
    from numpy.linalg import norm
    from scipy.constants import c

    h = zeros((N, N), dtype=complex)                        #Array of Green's scalars
    k0 = w0/c

    #Diagonal value is vacuum decay rate in units of gamma_0, must be set in hamiltonian "manually". See notes from meeting 7/3:
    diagval = 1

    if dimensionless:
        for i in range(N):
            for j in range(N):
                r = norm(rij[i,j])
                if i == j:
                    h[i, j] = diagval                             #Set diagonals to avoid divergence
                else:
                    h[i, j] = ((e**(1j * r)) / (4 * pi * r**3)) * (r**2 + 1j * r - 1)
    else:
        for i in range(N):
            for j in range(N):
                r = norm(rij[i,j])
                if i == j:
                    h[i, j] = diagval
                else:
                    h[i, j] = ((e**(1j * k0 * r)) / (4 * pi * k0**2 * r**3)) * (k0**2 * r**2 + 1j * k0 * r - 1)

    return h
                
