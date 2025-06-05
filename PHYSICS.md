# Effective Hamiltonian

A lattice of atoms situated closely together exercise strong interference effects in the electric field. For now, vacuum mediated E-field is chosen, and the following assumptions are made to scope in on the subspace that span the computational space:

 - The E-field is mediated much faster than the interaction, meaning retardation effects in relatively small lattices can be disregarded (in most cases),
 - Weak field coupling,
 - Atoms are identical and the same transition is chosen throughout, i.e. no detuning, which narrows the response bandwidth,
 - Tightly trapped atoms, i.e. stationary lattice sites.

Choosing the same transition in the atoms qubitifies the lattice, and the following notation for ground and excited states in the lattice is adopted (referred to as natural basis or atomic basis):

Ground state: $\ket{g} = \ket{g_0, g_1, \dots, g_{N-1}}$.

Singly excited state: $\ket{e_j} = \ket{g_0, \dots, g_{j-1}, e_j, g_{j+1}, \dots, g_{N-1}}$.

Doubly excited state: $\ket{e_j, e_k} = \ket{g_0, \dots, g_{j-1}, e_j, g_{j+1}, \dots, g_{k-1}, e_k, g_{k+1}, \dots, g_{N-1}}$.

Where N is the number of atoms in the lattice, and the e/g denotes if the indexed qubit is in the excited or ground state. 

In the dipole approximation, the output field is a superposition of the input field (e.g. a driving field) and the scattered field from the dipoles. This scattering off of dipoles can be described by a Green's tensor of classical field mediation in 3D, the dipole moment matrix element (i.e. transition strength) and a qubit deexcitation operator:

$$\hat{\bar{E}}_{out}(\bar{r}) = \hat{\bar{E}}_{in}(\bar{r}) + \mu_0 \omega_0^2 \sum_j \underline{G}(\bar{r}, \bar{r}_j, \omega_0) \bar{D}_j \hat{\sigma}_-^j$$.

Where $\mu_0$ is the vacuum permittivity, $\omega_0$ is the transition frequency, $\bar{r}_j$ is the j'th qubit position, $\bar{D}_j$ is the j'th dipole moment matrix elements of transition in 3D and $\hat{\sigma}_-^j$ is the j'th qubit deexcitation operator.

The i'th qubit can experience this outputted field, and the energy contribution of being excited from this field via dipole interaction gives the effective Hamiltonian: 

$$\hat{H}_{eff} = \sum_i -\hat{\bar{p}}_i \hat{\bar{E}}_{out}(\bar{r}_i) = - \mu_0 \omega_0^2 \sum_{i,j} \bar{D}_i^\dagger \underline{G}(\bar{r}_i, \bar{r}_j, \omega_0) \bar{D}_j \hat{\sigma}_+^i \hat{\sigma}_-^j$$.

In words, this is the energy contributions of the i'th qubit becoming excited from j'th qubit deexciting via dipole interactions with the i'th and j'th respective dipole moment matrix elements in the three spatial directions. The Hamiltonian is so called "effective", because it is not concerned with the complete system. The radiative modes, i.e. field modes that carry away the excitation from the lattice, are not included, which means every time-step, the probaility of exciting a such mode increases, and therefore the probability of still being excited decreases. This is reflected in the fact that the effective Hamiltonian is \textit{not} Hermitian. The Hamiltonian attains complex eigenvalues, where the imaginary value describes the eigenstate decay rate: 

$$\lambda_\xi = J_\xi - i \frac{\Gamma_\xi}{2}$$.

It is these eigenstates that can experience dramatically increased or decreased decay rates, known as super- and subradiant emission. The scope of interest is then: Can these states be identified and characterized? Can they be effectively adressed? Can the system be initialized in a such state, and can the system be manipulated to change these eigenvalues as desired? 

Working with non-Hermitian Hamiltonians prove difficult, because the spectral theorem is no longer valid. As discussed in Chapter 5 of [2], the quantum mechanics can be redeemed, and time-evolution of these eigenstates may still be carried out. 

In vacuum, the Green's tensor can be described ...

Considering also the effects of the qubits being excited on themselves, the full effective Hamiltonian is ...