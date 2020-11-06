# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import cirq
import scipy
import unicodedata
import numpy as np

# Define necessary Pauli operators (two-dimensional)
pauliX=np.array([[0,1],
                 [1,0]],
                dtype=complex)
pauliZ=np.array([[1,0],
                 [0,-1]],
                dtype=complex)

# Define dictionary of necessary Pauli operators (four-dimensional)
operator={
    "II": np.kron(np.identity(2),np.identity(2)),
    "IX": np.kron(np.identity(2),pauliX),
    "IZ": np.kron(np.identity(2),pauliZ),
    "XI": np.kron(pauliX,np.identity(2)),
    "XX": np.kron(pauliX,pauliX),
    "XZ": np.kron(pauliX,pauliZ),
    "ZI": np.kron(pauliZ,np.identity(2)),
    "ZX": np.kron(pauliZ,pauliX),
    "ZZ": np.kron(pauliZ,pauliZ)
    }

############################## State Preparation ##############################

''' 
Function that decomposes a tensor product of two unitary matrices into rotation 
gates, whose action is equivalent up to a global phase.
''' 
def tensorProductToRotationGates(unitaryMatrix0,unitaryMatrix1,qubit0,qubit1):
    
    # Compute the angles of each rotation from the entries of the matrix
    xa0=-np.angle(unitaryMatrix0[0][0])+np.angle(-unitaryMatrix0[0][1])
    xb0=2*np.arccos(np.abs(unitaryMatrix0[0][0]))
    xc0=-np.angle(unitaryMatrix0[0][0])+np.angle(unitaryMatrix0[1][0])
    xa1=-np.angle(unitaryMatrix1[0][0])+np.angle(-unitaryMatrix1[0][1])
    xb1=2*np.arccos(np.abs(unitaryMatrix1[0][0]))
    xc1=-np.angle(unitaryMatrix1[0][0])+np.angle(unitaryMatrix1[1][0])
    
    # Create the necessary rotation gates
    ra0=cirq.rz(xa0)
    rb0=cirq.ry(xb0)
    rc0=cirq.rz(xc0)
    ra1=cirq.rz(xa1)
    rb1=cirq.ry(xb1)
    rc1=cirq.rz(xc1)
    
    yield [ra0(qubit0),ra1(qubit1)]
    yield [rb0(qubit0),rb1(qubit1)]
    yield [rc0(qubit0),rc1(qubit1)]

''' 
Function that prepares a two qubit state passed as a 2x2 matrix, the 
element in position i,j being the coefficient of |ij>.
Returns a tuple with the operations to be performed first (ops), and the
single qubit unitary matrices that should be applied on the qubits afterwards
(unitary0 on qubit 0, unitary1 on qubit1).
''' 
def schmidtPreparationCircuit(stateMatrix,qubits):
    
    # Get the SVD decomposition
    unitary0, s, unitary1= np.linalg.svd(stateMatrix)
    
    # Create operations to get s[0]|00>+s[1]|11> from |00>
    rotation=cirq.ry(2*np.arccos(np.abs(s[0])))
    ops=[rotation(qubits[0]),cirq.CNOT(qubits[0],qubits[1])]
    
    return (ops,unitary0,unitary1)

'''  
Function that returns a circuit for the preparation of the two qubit state
specified by the input parameters.
''' 
def statePreparationCircuit(theta0,theta1,theta2,w0,w1,w2,qubits):
    
    # Get value of parametrized state coordinates
    alpha=np.cos(theta0/2)*np.cos(theta1/2)
    beta=np.cos(theta0/2)*np.sin(theta1/2)*np.exp(w1*1j)
    gamma=np.sin(theta0/2)*np.exp(w0*1j)*np.cos(theta2/2)
    delta=np.sin(theta0/2)*np.exp(w0*1j)*np.sin(theta2/2)*np.exp(w2*1j)
    
    matrix=np.array([[alpha,beta],
                    [gamma,delta]],
                    dtype=complex)

    ops,unitary0,unitary1=schmidtPreparationCircuit(matrix,qubits)
    yield ops
    yield tensorProductToRotationGates(unitary0,np.matrix.transpose(unitary1),\
                                qubits[0],qubits[1])

'''  
Function that returns a matrix to prepare a certain 2 qubit state from |00>.
Arguments: the six independent parameters that define the state.
Just for debugging.
 ''' 
def statePreparationMatrix(theta0,theta1,theta2,w0,w1,w2):
    # Get value of parametrized state coordinates
    alpha=np.cos(theta0/2)*np.cos(theta1/2)
    beta=np.cos(theta0/2)*np.sin(theta1/2)*np.exp(w1*1j)
    gamma=np.sin(theta0/2)*np.exp(w0*1j)*np.cos(theta2/2)
    delta=np.sin(theta0/2)*np.exp(w0*1j)*np.sin(theta2/2)*np.exp(w2*1j)
    
    # Create non-unitary matrix with the desired first column
    matrix=np.array([[alpha,1,0,0],
                     [beta,0,1,0],
                     [gamma,0,0,1],
                     [delta,0,0,0]],
                    dtype=complex)
    
    # Apply Gram-Schmidt process to create unitary matrix with
    #alpha, beta, gamma, delta as elements of the first column
    for i in range(1, len(matrix)):
        col_i = matrix[:, i]
        for j in range(0, i):
            col_j = matrix[:, j]
            t = col_i.dot(np.conj(col_j))
            col_i = col_i - t * col_j
        matrix[:, i] = col_i/np.linalg.norm(col_i)
        
    # Return matrix to prepare the state alpha, beta, gamma, delta from |00>
    return matrix

##################### Experimental Expectation Estimation #####################

#cirq.InsertStrategy.EARLIEST

''' 
Function that returns the expectation value of a given Pauli string.
Arguments: the string, the number of repetitions, the parameters defining 
the state in which to obtain the expectation value.
''' 

def measureExpectation(string,repetitions,theta0,theta1,theta2,w0,w1,w2):
    
    # Initialize qubits and circuit.
    n=2
    qubits=cirq.LineQubit.range(n)
    circuit=cirq.Circuit()
    
    # Append to the circuit the gates that prepare the state corresponding to
    #the received parameters.
    circuit.append(statePreparationCircuit(theta0,theta1,theta2,w0,w1,w2,qubits))
    
    # Append necessary rotations and measurements for each qubit.
    for i in range(n):
        op=string[i]
        
        # Rotate qubit i to the X basis if that's the desired measurement.
        if (op=="X"):
            circuit.append(cirq.H(qubits[i]))
            
        # Should include rotation to the Y basis if the Hamiltonian includes 
        #Pauli strings with Y matrix.
        #if (op=="Y")...
            
        # Measure qubit i in the computational basis, unless operator is I.
        if (op!="I"):
            circuit.append(cirq.measure(qubits[i],key=i))
            
    # Sample the desired number of repetitions from the circuit, unless
    #there are no measurements (II term).
    if (string!="II"):
        s=cirq.Simulator()
        results=s.run(circuit,repetitions=repetitions)
        
    # Calculate the expectation value of the Pauli string by averaging over  
    #all the repetitions.
    total=0
    for j in range(repetitions):
        meas=1
        for i in range(n):
            if (string[i]!="I"):
                meas=meas*(1-2*results.data[i][j])
        total+=meas
    expectationValue=total/repetitions
    return(expectationValue)

count=0

''' 
Function that returns the experimental energy expectation in a given state.
Arguments: the parameters defining the state in which to obtain the 
expectation value, a dictionary with the coefficient of each Pauli string
term in the Hamiltonian, the number of repetitions.
''' 
def experimentalExpectationEstimation(args,pauliStringCoefficient,repetitions):
    
    global count, maxfev
    if (count%(maxfev/10)==0):
        print(round(count/(maxfev/100)),"%",sep='')
    count=count+1
    
    # Get the parameters defining the state from the received arguments.
    theta0=args[0]
    theta1=args[1]
    theta2=args[2]
    w0=args[3]
    w1=args[4]
    w2=args[5]

    experimentalEnergyExpectation=0
    
    # Obtain experimental expectation value for each necessary Pauli string by
    #calling the measureExpectation function, and perform the necessary weighed
    #sum to obtain the energy expectation value.
    for pauliString in pauliStringCoefficient:
         
        expectationValue=measureExpectation(pauliString,repetitions,\
                                       theta0,theta1,theta2,w0,w1,w2)
        experimentalEnergyExpectation +=\
            pauliStringCoefficient[pauliString]*expectationValue

    return experimentalEnergyExpectation
    
##################### Theoretical Expectation Estimation #####################

''' 
Function that returns the theoretical energy expectation in a given state.
Arguments: the parameters defining the state in which to obtain the 
expectation value, a dictionary with the coefficient of each Pauli string
term in the Hamiltonian, the number of repetitions.
''' 
def theoreticalExpectationEstimation(args,pauliStringCoefficients):
     
    # Get the parameters defining the state from the received arguments.
    theta0=args[0]
    theta1=args[1]
    theta2=args[2]
    w0=args[3]
    w1=args[4]
    w2=args[5]

    theoreticalEnergyExpectation=0
    
    # Obtain theoretical expectation value for each necessary Pauli string by
    #matrix multiplication, and perform the necessary weighed sum to obtain 
    #the energy expectation value.
    for pauliString in pauliStringCoefficients:
    
        alpha=np.cos(theta0/2)*np.cos(theta1/2)
        beta=np.cos(theta0/2)*np.sin(theta1/2)*np.exp(w1*1j)
        gamma=np.sin(theta0/2)*np.exp(w0*1j)*np.cos(theta2/2)
        delta=np.sin(theta0/2)*np.exp(w0*1j)*np.sin(theta2/2)*np.exp(w2*1j)
        
        ket=np.array([alpha,beta,gamma,delta],dtype=complex)
        bra=np.conj(ket)
        
        expectationValue=np.real\
            (np.dot(bra,np.matmul(operator[pauliString],ket)))
        
        theoreticalEnergyExpectation+=\
            pauliStringCoefficients[pauliString]*expectationValue
    
    return theoreticalEnergyExpectation

############################# Complete Algorithm #############################

''' 
Function that uses the VQE algorithm to find the ground state and ground
energy of a Hamiltonian (2 qubits).
Arguments: the initial guess for the parameters, the Hamiltonian as a 
dictionary whose keys are Pauli strings and values the respective 
coefficients, the number of repetitions for the expectation estimation,
and a boolean flag 'simulate' that should be set to False to compute the
theoretical result, with no circuit simulations (for result checking).
''' 
def vqe(initialParameters,pauliStringCoefficients,repetitions=1000,\
        simulate=True):
    
    options={
        "disp": True,
        "maxfev": maxfev, # Maximum function evaluations
        "fatol": 0.0001, # Acceptable absolute error in f for convergence
        "xatol": 0.0001 # Acceptable absolute error in xopt for convergence
        }
    
    # Optimize the results from the CIRQ simulation
    if(simulate):
        optResults=scipy.optimize.minimize(experimentalExpectationEstimation,\
                                initialParameters,\
                                (pauliStringCoefficients,repetitions),\
                                    method='Nelder-Mead',options=options)
    
    # Optimize the results from theoretical computation
    else:
        optResults=scipy.optimize.minimize(theoreticalExpectationEstimation,\
                                initialParameters,(pauliStringCoefficients),\
                                    method='Nelder-Mead',\
                                    options=options)

    #print(optResults.message)
    # Print final parameters, obtained from optimization
    for i in range(6):
        
        if i<3:
            print(unicodedata.lookup("GREEK SMALL LETTER THETA"),i,": ",\
                  sep="",end="")
        else:
            print("w",i-3,": ",sep="",end="")
    
        # Bring parameters to [0,2pi] interval
        optimizedParameter=optResults.x[i]/np.pi
        while (optimizedParameter<0):
            optimizedParameter=optimizedParameter+2
        while (optimizedParameter>2):
            optimizedParameter=optimizedParameter-2
            
        print(optimizedParameter,unicodedata.lookup("GREEK SMALL LETTER PI")\
              ,sep="")
            
            
# Choose starting guess for the parameters.
theta0=0.2
theta1=2.1
theta2=1.3
w0=4.3
w1=5.8
w2=1.2
    
# Define the coefficient of each Pauli string term in the Hamiltonian   
# Here R=0.90
pauliStringCoefficients90={
    "II": -3.8505,
    "IX": -0.2288,
    "IZ": -1.0466,
    "XI": -0.2288,
    "XX": +0.2613,
    "XZ": +0.2288,
    "ZI": -1.0466,
    "ZX": +0.2288,
    "ZZ": +0.2356
    }

maxfev=1000 # Maximum function evaluations (for the optimization)

#vqe([theta0,theta1,theta2,w0,w1,w2],pauliStringCoefficients90,simulate=False)

#################################### Tests ####################################

# Generates random paramenters, prepares the corresponding state, 
#and returns a bool that indicates if the state preparation is fine.
def checkRandomStatePreparation():
    theta0=np.pi*np.random.rand()
    theta1=np.pi*np.random.rand()
    theta2=np.pi*np.random.rand()
    w0=2*np.pi*np.random.rand()
    w1=2*np.pi*np.random.rand()
    w2=2*np.pi*np.random.rand()
    
    matrix=statePreparationMatrix(theta0,theta1,theta2,w0,w1,w2)
    qubits=cirq.LineQubit.range(2)
    circuit=cirq.Circuit()
    circuit.append(cirq.MatrixGate(matrix).on(qubits[0],qubits[1]))
    s=cirq.Simulator()
    results=s.simulate(circuit)
    
    circuit2=cirq.Circuit()
    qubits2=cirq.LineQubit.range(2)
    circuit2.append(statePreparationCircuit(theta0,theta1,theta2,w0,w1,w2,qubits2))
    s=cirq.Simulator()
    results2=s.simulate(circuit2)
    
    return(cirq.equal_up_to_global_phase(results.state_vector(),\
                                         results2.state_vector()))
        
def checkStatePreparation(numberOfTests):
    k=0
    for i in range(numberOfTests):
        k=k+checkRandomStatePreparation()
    return k

numberOfTests=100
k=checkStatePreparation(numberOfTests)
print("For",k,"out of",numberOfTests,\
      "random tests, the state preparations were identical up to a global phase.")

