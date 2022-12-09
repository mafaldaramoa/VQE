# Ansätze for noisy variational quantum eigensolvers
Mafalda Ramôa

This repository contains the code I developed for my thesis, submitted for the degree of Master of Science in the Engineering Physics course at University of Minho (Physics of Information branch).

[You can access the thesis on arXiv.](https://arxiv.org/abs/2212.04323)

## Files

### CIRQ_AdaptVQE

Implementation of Adapt VQE ([1][2]), using the CIRQ simulator for the evaluation of both the energy and the gradient, and OpenFermion, PySCF and their plug-in for the chemistry part.
Includes:
- A noise free, matrix algebra based version of the algorithm.
- The possibility of removing terms or growing the ansatz conservatively (multiple optimization attempts per iteration).
- Multiple pool options (the ones used in Qubit Adapt [1] and Fermionic Adapt [2], and others).
- Functions to estimate the required numer of CNOTs.
- Functions to examine the Slater determinants in the state.
- Plots for the evolution of the energy, comparing runs, etc.

### CIRQ_AdaptVQE_Basic
Simpler, function-based version of Adapt VQE using the CIRQ simulator  for the evaluation of both the energy and the gradient, and OpenFermion, PySCF and their plug-in for the chemistry part.
Includes:
- A noise free, matrix algebra based version of the algorithm.
- Multiple pool options (the ones used in Fermionic Adapt and Qubit Adapt, and others).
- Energy calculation functions, with and without sampling noise, and tests to check their performance, as well as the impact of trotterization and sampling noise on the energy evaluation.
- Gradient calculation functions, with and without sampling noise, and tests to check their performance, as well as the impact of the number of shots on the precision.
- Basic energy plots.

### CIRQ_Openfermion 
Exploring simple functionalities of Openfermion along with CIRQ.

### CIRQ_UCCSD_VQE
Implementation of VQE with the unitary coupled cluster singles and doubles ansatz, using the CIRQ simulator for the energy expectation, and Openfermion, PySCF and their plug-in for the chemistry part.
Includes:
- A noise free, matrix algebra based version of the algorithm.
- Analysis of the number of operators in the UCCSD ansatz.
- Tests on the expectation estimation, with and without trotterization, and with or without sampling noise.
- Multiple methods for trotterizing operators, and tests on their application to Hamiltonians and the UCCSD operator.
- Analysis of the effect of the starting point, trotterization, sampling noise, and choice of optimizer.

### CIRQ_VQE 
Implementation of the original VQE algorithm [3] using the CIRQ simulator for the energy expectation. The algorithm is applied to HeH+, using the Hamiltonian from that paper, and the same ansatz (spanning the whole 2-qubit Hilbert space).
Includes:
- A noise free, matrix algebra based version of the algorithm.
- Functions to create a circuit to prepare an arbitrary two qubit state from its parameterization, based on Schmidt decomposition.
- Plot comparing the median and average of the energy.
- Plot with the median of the energy, with interquartile ranges as error bars.
- Tests on the energy expectation and state preparation.

### CIRQ_VQE_Hyperparameters
Attempts to minimize the hyperparameters for the Nelder Mead optimization: tolerance (fatol, xatol) and the initial simplex size (through a single parameter delta). 
Current values were result of using a cost function penalizing overlaps with the ground state under 90%, and the number of failed attempts until converging.

### Qiskit_AdaptVQE
Implementation of the Adapt VQE algorithm, using Qiskit for evaluation of the energy, and OpenFermion, PySCF and their plug-in for the chemistry part.
- A noise free, matrix algebra based version of the algorithm..
- Noise models for thermal relaxation, SPAM noise,... 
- Qubit Adapt [1] and Fermionic Adapt [2] pools.
- Plot for the evolution of the energy
- Trotterization of OpenFermion's QubitOperator into Qiskit circuits.
- Transformation of OpenFermion's Hamiltonian into Qiskit Aqua observables.


### Qiskit_UCCSD_VQE
Implementation of VQE with the unitary coupled cluster singles and doubles ansatz, using Qiskit for the energy expectation, and PySCF and its Qiskit Driver for the chemistry part.
Includes:
- A noise free, matrix algebra based version of the algorithm..
- Noise models for thermal relaxation, SPAM noise,... 
- Numerical results.
- A plot of the evolution of the energy along the optimization.

### Qiskit_VQE
Implementation of the original VQE algorithm [3] using Qiskit for the energy expectation. The algorithm is applied to HeH+, using the Hamiltonian from that paper, and the same ansatz (spanning the whole 2-qubit Hilbert space).
Includes:
- A noise free, matrix algebra based version of the algorithm.
- Plot of the evolution of the energy and overlap with the ground state along the optimization.

## References
[1]  Ho  Lun  Tang,  V.  O.  Shkolnikov, George S.  Barron, Harper R.  Grimsley, Nicholas J. Mayhall, Edwin Barnes, and Sophia E. Economou. qubit-ADAPT-VQE:  An adaptive algorithm for constructing hardware-efficient ansatze on a quantum processor. Preprint arXiv arXiv:1911.10205 [quant-ph].

[2] Harper R. Grimsley, Sophia E. Economou, Edwin Barnes, and Nicholas J. Mayhall. An adaptive variational algorithm for exact molecular simulations on a quantum computer. Nature Communications, 10(1), 3007 (2019).

[3] Alberto Peruzzo, Jarrod McClean, Peter Shadbolt, Man-Hong Yung, Xiao-Qi Zhou, Peter J.  Love, Al ́an Aspuru-Guzik, and Jeremy L.  O’Brien.   A variational eigenvalue solver on a photonic quantum processor. Nature Communications 5 (1), 4213 (2014).
