# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import cirq
import scipy
import unicodedata
import numpy as np
import statistics
import matplotlib.pyplot as plt
import xlrd

# Define necessary Pauli operators (two-dimensional) as matrices
pauliX=np.array([[0,1],
                 [1,0]],
                dtype=complex)
pauliZ=np.array([[1,0],
                 [0,-1]],
                dtype=complex)

# Define dictionary of necessary Pauli operators (four-dimensional) as matrices
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

# Example Hamiltonian, for atomic distance R=90pm
hamiltonian90={
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

# Read the data table with the coefficients of Pauli strings for the 
#Hamiltonians of several atomic distances for He-H+
#(article by Peruzzo et al)
directory=r'C:\Users\mafal\OneDrive\Computador\Universidade\Tese(VQE)\DataPeruzzo.xlsx'
workbook = xlrd.open_workbook(directory)
    
worksheet = workbook.sheet_by_index(0)

# Store in a list the names of the Pauli strings that appear in the first row
pauliStrings = []

for row in range(worksheet.ncols):
    pauliStrings.append( worksheet.cell_value(0,row) )
    
# Create a list of Hamiltonians as dictionaries from the data, and a list of
#the atomic distances
listOfHamiltonians=[]
listOfRadii=[]

for row in range(1, worksheet.nrows):
    hamiltonian={}
    
    # Append to the list of atomic separations the value in picometers
    listOfRadii.append(100*float(worksheet.cell_value(row,0))) 
    
    # Associate to each Pauli string (key) in the Hamiltonian (dictionary) 
    #its coefficient (value)
    for col in range(1,worksheet.ncols):
        hamiltonian[pauliStrings[col]]=float(worksheet.cell_value(row,col))
        
    # Append the dictionary representing the Hamiltonian to the list of 
    #Hamiltonians
    listOfHamiltonians.append(hamiltonian)

############################## State Preparation ##############################

def fromParametersToCoordinates(stateParameters):
    '''
    Arguments: a list containing 6 parameters defining a normalized 2 qubit state, 
    with the DoF corresponding to the global phase removed.
    Returns: a list containing the 4 complex coordinates of the state in the 
    computational basis.
    '''
    
    # Get the parameters from the input list
    theta0=stateParameters[0]
    theta1=stateParameters[1]
    theta2=stateParameters[2]
    w0=stateParameters[3]
    w1=stateParameters[4]
    w2=stateParameters[5]
    
    # Calculate complex coordinates from parameters
    alpha=np.cos(theta0/2)*np.cos(theta1/2)
    beta=np.cos(theta0/2)*np.sin(theta1/2)*np.exp(w1*1j)
    gamma=np.sin(theta0/2)*np.exp(w0*1j)*np.cos(theta2/2)
    delta=np.sin(theta0/2)*np.exp(w0*1j)*np.sin(theta2/2)*np.exp(w2*1j)
    
    return [alpha,beta,gamma,delta]

def tensorProductToRotationGates(unitaryMatrix0,unitaryMatrix1,qubit0,qubit1,\
                                 keepPhase):
    ''' 
    Arguments: two unitary matrices, two qubits, a boolean keepPhase indicating
    whether global phase should be kept.
    Returns: a circuit composed of single qubit rotation gates, that is 
    equivalent (up to a global phase) to the action of the tensor product of 
    the two matrices. 
    ''' 
    
    # Compute the angles of each rotation from the entries of the matrix
    rz1angle_q0=-np.angle(unitaryMatrix0[0][0])+np.angle(-unitaryMatrix0[0][1])
    ryangle_q0=2*np.arccos(np.abs(unitaryMatrix0[0][0]))
    rz2angle_q0=-np.angle(unitaryMatrix0[0][0])+np.angle(unitaryMatrix0[1][0])
    
    rz1angle_q1=-np.angle(unitaryMatrix1[0][0])+np.angle(-unitaryMatrix1[0][1])
    ryangle_q1=2*np.arccos(np.abs(unitaryMatrix1[0][0]))
    rz2angle_q1=-np.angle(unitaryMatrix1[0][0])+np.angle(unitaryMatrix1[1][0])
    
    if keepPhase:
        gp0=np.angle(unitaryMatrix0[0][0])+(rz1angle_q0+rz2angle_q0)/2
        gp1=np.angle(unitaryMatrix1[0][0])+(rz1angle_q1+rz2angle_q1)/2
        gp=gp0+gp1
    
    # Create the necessary rotation gates
    r1_q0=cirq.rz(rz1angle_q0+rz1angle_q1)
    r2_q0=cirq.ry(ryangle_q0)
    r3_q0=cirq.rz(rz2angle_q0)
    
    r1_q1=cirq.ry(ryangle_q1)
    r2_q1=cirq.rz(rz2angle_q1)
    
    '''
    Reconstruct matrix from the calculated angles.
    Just for debbugging.
    if(keepPhase):
        
        # Reconstruct matrix acting on qubit 0
        u00=np.cos(ryangle_q0/2)*np.exp(-1j*((rz1angle_q0+rz2angle_q0)/2-gp0))
        u10=np.sin(ryangle_q0/2)*np.exp(-1j*((rz1angle_q0-rz2angle_q0)/2-gp0))
        u01=-np.sin(ryangle_q0/2)*np.exp(-1j*((-rz1angle_q0+rz2angle_q0)/2-gp0))
        u11=np.cos(ryangle_q0/2)*np.exp(-1j*((-rz1angle_q0-rz2angle_q0)/2-gp0))
        print(u00,u01,u10,u11)
        print(unitaryMatrix0)
        
        # Reconstruct matrix acting on qubit 1
        v00=np.cos(ryangle_q1/2)*np.exp(-1j*((rz1angle_q1+rz2angle_q1)/2-gp1))
        v10=np.sin(ryangle_q1/2)*np.exp(-1j*((rz1angle_q1-rz2angle_q1)/2-gp1))
        v01=-np.sin(ryangle_q1/2)*np.exp(-1j*((-rz1angle_q1+rz2angle_q1)/2-gp1))
        v11=np.cos(ryangle_q1/2)*np.exp(-1j*((-rz1angle_q1-rz2angle_q1)/2-gp1))
        print(v00,v01,v10,v11)
        print(unitaryMatrix1)
    '''
    
    yield [r1_q0(qubit0)]
    yield [r2_q0(qubit0),r1_q1(qubit1)]
    yield [r3_q0(qubit0),r2_q1(qubit1)]
    if keepPhase:
        yield cirq.MatrixGate(np.identity(2)*np.exp(1j*gp)).on(qubit0)
    
def schmidtDecompositionCircuit(stateMatrix,qubits):
    ''' 
    Arguments: a 2x2 matrix representing a two qubit state, the element in 
    position i,j being the coefficient of |ij>, and a pair of qubits to be 
    prepared in that state.
    Returns: a recipe to prepare the state, as a tuple. "gates" is the list of
    gates to be applied first to create the necessary amount of entanglement, 
    and unitary0, unitary1 are the single qubit unitary matrices representing 
    the operations that should be applied  on the qubits (in position 0 and 1 
    respectively) afterwards.
    ''' 
    
    # Get the SVD decomposition
    unitary0, s, unitary1= np.linalg.svd(stateMatrix)
    
    # Create operations to get s[0]|00>+s[1]|11> from |00>
    rotation=cirq.ry(2*np.arccos(np.abs(s[0])))
    gates=[rotation(qubits[0]),cirq.CNOT(qubits[0],qubits[1])]
    
    return (gates,unitary0,np.matrix.transpose(unitary1))

def statePreparationGates(stateParameters,qubits,keepPhase=False):
    '''  
    Returns: gates to be applied on two qubits to prepare a certain state from 
    |00>.
    Arguments: parameters defining  the two qubit state, and the qubits to be 
    prepared on this state.
    ''' 
    
    # Get the state coordinates
    alpha,beta,gamma,delta=fromParametersToCoordinates(stateParameters)
    
    # Write the state coordinates as a 2x2 matrix
    matrix=np.array([[alpha,beta],
                    [gamma,delta]],
                    dtype=complex)

    # Get the gates to be applied first, and the unitary single qubit matrices
    #to be applied second to prepare the state
    gates,unitary0,unitary1=schmidtDecompositionCircuit(matrix,qubits)
    
    yield gates
    yield tensorProductToRotationGates(unitary0,unitary1,qubits[0],qubits[1],\
                                       keepPhase)

##################### Experimental Expectation Estimation #####################


def measureExpectation(pauliString,repetitions,stateParameters):
    ''' 
    Returns: the expectation value of a Pauli string, obtained by using the
    CIRQ simulator.
    Arguments: the Pauli string to be measured, the number of repetitions to be 
    performed, the list of parameters defining the state in which to obtain the 
    expectation value.
    ''' 
    
    # Initialize qubits and circuit.
    n=2
    qubits=cirq.LineQubit.range(n)
    circuit=cirq.Circuit()
    
    # Append to the circuit the gates that prepare the state corresponding to
    #the received parameters.
    circuit.append(statePreparationGates(stateParameters,qubits))
    
    # Append necessary rotations and measurements for each qubit.
    for i in range(n):
        op=pauliString[i]
        
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
    #there are no measurements (identity term).
    if (pauliString!="II"):
        s=cirq.Simulator()
        results=s.run(circuit,repetitions=repetitions)
        
    # Calculate the expectation value of the Pauli string by averaging over  
    #all the repetitions.
    
    total=0
    
    for j in range(repetitions):
        meas=1
        for i in range(n):
            if (pauliString[i]!="I"):
                meas=meas*(1-2*results.data[i][j])
        total+=meas
        
    expectationValue=total/repetitions
    
    return(expectationValue)

count=0
optimizing=False

def experimentalExpectationEstimation(stateParameters,pauliStringCoefficients\
                                      ,repetitions):
    ''' 
    Returns: the experimental energy expectation in a given state.
    Arguments: a list of 6 parameters defining the state in which to obtain the 
    expectation value, a dictionary with the coefficients of each Pauli string
    term in the Hamiltonian, the number of repetitions to be used to calculate 
    the expectation value of each string.
    ''' 
    
    # Print the percentage of the maximum number of function evaluations that 
    #has been used so far in the classical optimization, if it's a multiple of
    #10% (to inform on the stage of the optimization process).
    global count, optimizing
    if(optimizing):
        if (count%(300/10)==0):
            print(round(count/(300/100)),"%",sep='')
        count=count+1

    experimentalEnergyExpectation=0
    
    # Obtain experimental expectation value for each necessary Pauli string by
    #calling the measureExpectation function, and perform the necessary weighed
    #sum to obtain the energy expectation value.
    for pauliString in pauliStringCoefficients:
         
        expectationValue=measureExpectation(pauliString,repetitions,\
                                       stateParameters)
        experimentalEnergyExpectation +=\
            pauliStringCoefficients[pauliString]*expectationValue

    return experimentalEnergyExpectation
    
##################### Theoretical Expectation Estimation #####################

def theoreticalExpectationEstimation(stateParameters,pauliStringCoefficients):
    ''' 
    Returns: the theoretical energy expectation in a given state.
    Arguments: a list with the 6 parameters defining the state in which to 
    obtain the expectation value, a dictionary with the coefficients of each
    Pauli string term in the Hamiltonian.
    ''' 

    theoreticalEnergyExpectation=0
    
    # Obtain the theoretical expectation value for each Pauli string in the
    #Hamiltonian by matrix multiplication, and perform the necessary weighed
    #sum to obtain the energy expectation value.
    for pauliString in pauliStringCoefficients:
    
        alpha,beta,gamma,delta=fromParametersToCoordinates(stateParameters)
        
        ket=np.array([alpha,beta,gamma,delta],dtype=complex)
        bra=np.conj(ket)
        
        expectationValue=np.real\
            (np.dot(bra,np.matmul(operator[pauliString],ket)))
        
        theoreticalEnergyExpectation+=\
            pauliStringCoefficients[pauliString]*expectationValue
            
    return theoreticalEnergyExpectation

############################# Complete Algorithm #############################

def vqe(initialParameters,pauliStringCoefficients,repetitions=1000,\
        simulate=True):
    ''' 
    Returns: an OptimizeResult object consisting of the result of attempting to
    minimize the expectation value of the energy by using the VQE algorithm.
    Arguments: a list with the initial guess for the 6 parameters, the 
    Hamiltonian as a dictionary whose keys are Pauli strings and values the 
    respective coefficients, the number of repetitions for the expectation 
    estimation, and a boolean flag 'simulate' that should be set to False to 
    compute the theoretical result, with no circuit simulations 
    (for result checking).
    ''' 
    global optimizing 
    
    # Choose maximum number of function evaluations for the optimization
    # A lower number seems to work better when running the CIRQ simulator
    if(simulate):
        maxfev=300
    else:
        maxfev=1000
    
    # Select the options for the optimization
    options={
        "disp": True,
        "maxfev": maxfev, # Maximum function evaluations
        "fatol": 0.05, # Acceptable absolute error in f for convergence
        "xatol": 0.1, # Acceptable absolute error in xopt for convergence
        "adaptive": True
        }
    
    optimizing=True
    
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
                                    method='Nelder-Mead',
                                    options=options)
    optimizing=False
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
   
    return optResults

############################## Testing Functions ##############################

def statePreparationMatrix(stateParameters):
    '''  
    Returns: a matrix to prepare a certain 2 qubit state from |00>.
    Arguments: a list containing the six parameters that define the state.
    Just for debugging.
    ''' 
    
    # Get value of parametrized state coordinates
    alpha,beta,gamma,delta=fromParametersToCoordinates(stateParameters)
    
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
        
    return matrix

def testExpectationEstimation(repetitions,hamiltonian):
    '''
    Returns: a tuple with the error and the standard deviation of the quantum 
    #expectation estimation using the CIRQ simulator (taking the average over
    the results for 10 random states, 100 QEE runs in each).
    Arguments: the number of repetitions used in QEE.
    '''
        
    stdList=[]
    errorList=[]
    
    # Calculate standard deviation and error for 10 random states
    for _ in range(10):
        
        # Generate list of 6 parameters defining a random 2 qubit state 
        theta0=np.pi*np.random.rand()
        theta1=np.pi*np.random.rand()
        theta2=np.pi*np.random.rand()
        w0=2*np.pi*np.random.rand()
        w1=2*np.pi*np.random.rand()
        w2=2*np.pi*np.random.rand()
        stateParameters=[theta0,theta1,theta2,w0,w1,w2]
        
        # Obtain the theoretical value for the energy
        theoreticalEnergyExpectation=theoreticalExpectationEstimation\
            (stateParameters,hamiltonian)
        
        expectationList=[]
        
        # Obtain the experimental value for the energy in 100 QEE runs
        for _ in range(100):
                    
            # Test the experimental energy expectation against the theoretical 
            #value
            experimentalEnergyExpectation=experimentalExpectationEstimation\
                (stateParameters,hamiltonian,repetitions)
            
            expectationList.append(experimentalEnergyExpectation)
            
        # Append to the error list the average difference between 
        #experimental and theoretical energy over the 100 runs
        error=np.mean(np.abs(expectationList-theoreticalEnergyExpectation))
        errorList.append(error)
        
        # Append to the standard deviation list the standard deviation of the
        #100 runs
        stdList.append(np.std(expectationList))
    
    # Return the error and standard deviation estimated by taking the average
    #for the 10 random states
    avError=np.mean(errorList)
    avSTD=np.mean(stdList)
    
    return (avError,avSTD)

def pauliStringSTD(pauliString,repetitions):
    '''
    Returns: the standard deviation associated with the expectation value of a 
    given Pauli string (taking the average over the results for 10 random 
    states, 100 QEE runs in each).
    Arguments: the Pauli string in question, and the number of repetitions for 
    QEE.
    '''

    stdList=[]
    
    # Calculate standard deviation for 10 random states
    for _ in range(10):
        
        # Generate list of 6 parameters defining a random 2 qubit state 
        theta0=np.pi*np.random.rand()
        theta1=np.pi*np.random.rand()
        theta2=np.pi*np.random.rand()
        w0=2*np.pi*np.random.rand()
        w1=2*np.pi*np.random.rand()
        w2=2*np.pi*np.random.rand()
        stateParameters=[theta0,theta1,theta2,w0,w1,w2]
        
        # Get the expectation estimation for the Pauli string in the same state
        #in 100 different runs
        expectationList=[]
        for _ in range(100):
                    
            experimentalExpectation=measureExpectation\
                (pauliString,repetitions,stateParameters)
            
            expectationList.append(experimentalExpectation)
        
        # Add the standard deviation resulting from this state to the list 
        stdList.append(np.std(expectationList))
        
    # Return the average of the standard deviation over the 10 different 
    #states
    return np.mean(stdList)

def dictOfSTD(repetitions,listOfPauliStrings):
    '''
    Returns: a dictionary of the standard deviation associated with the 
    expectation value of multiple Pauli strings, obtained via the CIRQ 
    simulator.
    Arguments: the number of repetitions to be used in obtaining the 
    expectation value, a list containing the Pauli strings whose standard 
    deviation will be calculated.
    '''
    
    STDDictionary={}
    
    for pauliString in listOfPauliStrings:
        STDDictionary[pauliString]=pauliStringSTD(pauliString,repetitions)
        
    return STDDictionary
    

def errorPropagation(stdDictionary,hamiltonian):
    '''
    Returns: the error in the energy expecation value obtained via the CIRQ 
    simulator, resulting from propagating the error associated with the 
    multiple Pauli strings in a given Hamiltonian.
    Arguments: a dictionary containing Pauli strings as keys (at least the ones
    included in the Hamiltonian) and the respective standard deviation as 
    values, and the Hamiltonian as a dictionary of Pauli strings and respective
    coefficients.
    '''
    
    error=0
    
    # Add to the error the variance of each term
    for pauliString in hamiltonian:
        error=error+pow(stdDictionary[pauliString]*hamiltonian[pauliString],2)
            
    # Take the square root to obtain the standard deviation
    error=np.sqrt(error)
    
    return error

def stateFromGates(gates,repetitions,qubits):
    '''
    Returns: a matrix with the element in row i, column j being the probability
    of obtaining state |ij> at the end of a given circuit.
    Arguments: the gates the circuit is comprised of, the number of repetitions 
    to be used to obtain the probabilities, the two qubits the circuit acts on.
    '''
    
    # Create a circuit from the provided gates
    circuit=cirq.Circuit(gates)
    
    # Add measurements and run the circuit
    circuit.append(cirq.measure(qubits[0],key=0))
    circuit.append(cirq.measure(qubits[1],key=1))
    s=cirq.Simulator()
    results=s.run(circuit,repetitions=repetitions)
    
    # Calculate the probability of occurence of each computational basis state
    total=np.array([[0,0],[0,0]])
    
    for j in range(repetitions):
        total[results.data[0][j]][results.data[1][j]]=\
            total[results.data[0][j]][results.data[1][j]]+1
            
    probabilityMatrix=total/repetitions
    
    # Return the probabilities organized as a matrix
    return(probabilityMatrix)
    
def testPreparationCircuit(numberOfTests,repetitions):
    '''
    Returns: the average error in probability amplitudes obtained by running
    on CIRQ the circuit created by statePreparationGates to prepare a state, 
    as compared to the actual probabilities, calculated from the wavefunction.
    Arguments: the number of tests over which the average error should be 
    calculated, and the number of repetitions to be used in obtaining the
    experimental probabilities.
    '''
    
    avError=0
    
    # Calculate the average error in probabilities over the desired number of 
    #tests
    for _ in range(numberOfTests):
        
        # Generate 6 random parameters defining a state
        theta0=np.pi*np.random.rand()
        theta1=np.pi*np.random.rand()
        theta2=np.pi*np.random.rand()
        w0=2*np.pi*np.random.rand()
        w1=2*np.pi*np.random.rand()
        w2=2*np.pi*np.random.rand()
        stateParameters=[theta0,theta1,theta2,w0,w1,w2]
    
        # Calculate the theoretical probability matrix (element i,j being 
        #the probability of state |ij>)
        alpha,beta,gamma,delta=fromParametersToCoordinates(stateParameters)
        
        p00=pow(np.abs(alpha),2)
        p01=pow(np.abs(beta),2)
        p10=pow(np.abs(gamma),2)
        p11=pow(np.abs(delta),2)
        
        theoreticalProbabilityMatrix=[[p00,p01],\
                                      [p10,p11]]
        
        # Obtain the experimental probability matrix
        qubits=cirq.LineQubit.range(2)
        experimentalProbabilityMatrix=stateFromGates(statePreparationGates\
                         (stateParameters,qubits),repetitions,qubits)
            
        absoluteErrorMatrix=np.abs(theoreticalProbabilityMatrix\
                                     -experimentalProbabilityMatrix)
            
        avError=avError+np.mean(absoluteErrorMatrix)
        
        return avError/numberOfTests
    
#################################### Tests ####################################

# Set flags to True to execute tests
testStatePreparation = False # State preparation
testQEE = False # Expectation estimation

# Choose the number of random tests to be performed on the state preparation
numberOfTests = 1000

# Choose the number of shots to be used in the state preparation tests
repetitions = 100

# Choose number of shots to be used in the expectation estimation tests
repetitionsQEE = 100

# Choose hamiltonian to be used in the expectation estimation tests
hamiltonian=hamiltonian90


if testStatePreparation:
    # Test the state preparation circuit

    avError=testPreparationCircuit(numberOfTests,repetitions)

    print("***Testing the state preparation circuit***")
    print("\nBy getting the statistics:")
    print("Completed",numberOfTests,"random tests, with",repetitions,\
          "shots.","\nAverage error in outcome probabilities:",avError)
        
    theta0=np.pi*np.random.rand()
    theta1=np.pi*np.random.rand()
    theta2=np.pi*np.random.rand()
    w0=2*np.pi*np.random.rand()
    w1=2*np.pi*np.random.rand()
    w2=2*np.pi*np.random.rand()
    
    stateParameters=[theta0,theta1,theta2,w0,w1,w2]
    
    # Calculate the theoretical coordinates
    coordinates=fromParametersToCoordinates(stateParameters)
    
    roundedCoordinates=np.ndarray(4,dtype=complex)
    for i,coordinates in enumerate(coordinates):
        roundedCoordinates[i]=(np.round(coordinates,7))
        
    qubits=cirq.LineQubit.range(2)
    circuit=cirq.Circuit()
    circuit.append(statePreparationGates(stateParameters,qubits,keepPhase=True))
    simulator=cirq.Simulator()
    s=simulator.simulate(circuit)
    
    print("\nBy simulating the full wavefunction:")
    print("Actual state coordinates:")
    print(roundedCoordinates)
    print("State resulting from the preparation circuit:")
    print(s.final_state_vector)

if testQEE:
    # Test the expectation estimation
    
    # By testing the full QEE function 
    error,std=testExpectationEstimation(repetitionsQEE,hamiltonian)
    
    print("***Testing the expectation estimation for",repetitions,"shots***")
    print("By testing the full function:")
    print("Obtained error",error,"and standard deviation",std)
    
    # By propagating the error from the individual terms
    listOfPauliStrings=hamiltonian.keys()
    dictOfSTD=dictOfSTD(repetitionsQEE,listOfPauliStrings)
    print("By propagating the error:")
    print("Obtained error",errorPropagation(dictOfSTD,hamiltonian))
    
    # Propagate error for all points in the bond dissociation curve
    errorList=[]
    
    # Calculate error for all the atomic distances
    for hamiltonian in listOfHamiltonians:
        errorList.append(errorPropagation(dictOfSTD,hamiltonian))
        
    # Plot the results
    plt.title("Error in the Expectation Estimation ("+str(repetitions)+" shots)")
    plt.xlabel("Atomic distance (pm)")
    plt.ylabel("Error (MJ mol$^{-1}$)")
    plt.plot(listOfRadii,errorList)

################## Functions for the Bond Dissociation Graph ##################

def calculateOverlap(stateCoordinates1,stateCoordinates2):
    '''
    Returns: the overlap between two states.
    Arguments: the coordinates defining each state in the computational basis.
    '''
    
    # Calculate the absolute value of the inner product of the two states
    overlap=np.abs(np.dot(np.conj(stateCoordinates1),stateCoordinates2))
    
    return overlap

def groundStatesFromDiagonalization(listOfHamiltonians):
    '''
    Returns: a tuple with a list of the ground energies and a list of the 
    ground states corresponding to a list of hamiltonians, obtained by exact 
    diagonalization.
    Arguments: the list of hamiltonians as dictionaries.
    '''
    
    exactGroundEnergies=[]
    exactGroundStates=[]
    
    for hamiltonian in listOfHamiltonians:
        hamiltonianMatrix=np.zeros((4,4),dtype=complex)
        
        for pauliString in hamiltonian:
            hamiltonianMatrix+=hamiltonian[pauliString]*operator[pauliString]
            
        w,v = np.linalg.eig(hamiltonianMatrix)
        
        exactGroundEnergies.append((np.amin(w)))
        exactGroundStates.append(v[np.argmin(w)])
        
    return(exactGroundEnergies,exactGroundStates)
            
def groundStatesFromTheoreticalOptimization(listOfHamiltonians,runs):
    '''
    Returns: a tuple with a list of the ground energies and a list of the 
    ground states corresponding to a list of hamiltonians, obtained by 
    optimizing the value coming from exact calculation of the energy.
    Arguments: the list of hamiltonians as dictionaries, and the number of runs
    over which to take the median (should be an odd number).
    '''
    
    groundEnergies=[]
    groundStates=[]
    
    for hamiltonian in listOfHamiltonians:
            
        optimizationRuns=[]
        stateArray=[]
        
        for i in range(runs):
            
            success=False
            while not success:
                
                theta0=np.pi*np.random.rand()
                theta1=np.pi*np.random.rand()
                theta2=np.pi*np.random.rand()
                w0=2*np.pi*np.random.rand()
                w1=2*np.pi*np.random.rand()
                w2=2*np.pi*np.random.rand()
                
                results=vqe([theta0,theta1,theta2,w0,w1,w2],hamiltonian,\
                            simulate=False)
                success=results.success
                
            groundEnergy=results.fun
            optimizationRuns.append(groundEnergy)
            stateArray.append(results.x)
            
        groundEnergy=statistics.median(optimizationRuns)
        groundEnergies.append(groundEnergy)
        groundState=fromParametersToCoordinates\
            (stateArray[optimizationRuns.index(groundEnergy)])
        groundStates.append(groundState)
        
    return(groundEnergies,groundStates)

def groundStatesFromVQE(listOfHamiltonians,repetitions,runs):
    '''
    Returns: a tuple with a list of the ground energies and a list of the 
    ground states corresponding to a list of hamiltonians, obtained by 
    optimizing the value coming from QEE on the CIRQ simulator.
    Arguments: the list of hamiltonians as dictionaries, the number of
    repetitions to be used in QEE, and the number of runs over which to take
    the median (should be an odd number).
    '''
    
    global count
     
    groundEnergies=[]
    groundStates=[]
    
    for hamiltonian in listOfHamiltonians:
            
        optimizationRuns=[]
        stateArray=[]
        
        for i in range(runs):
            
            success=False
            while not success:
                
                theta0=np.pi*np.random.rand()
                theta1=np.pi*np.random.rand()
                theta2=np.pi*np.random.rand()
                w0=2*np.pi*np.random.rand()
                w1=2*np.pi*np.random.rand()
                w2=2*np.pi*np.random.rand()
                
                count=0
                results=vqe([theta0,theta1,theta2,w0,w1,w2],hamiltonian,\
                        repetitions=repetitions)
                success=results.success
                
            groundEnergy=results.fun
            optimizationRuns.append(groundEnergy)
            stateArray.append(results.x)
            
        groundEnergy=statistics.median(optimizationRuns)
        groundEnergies.append(groundEnergy)
        groundState=fromParametersToCoordinates\
            (stateArray[optimizationRuns.index(groundEnergy)])
        groundStates.append(groundState)
        
    return(groundEnergies,groundStates)

########################### Bond Dissociation Graph ###########################

# Set flag to True to generate the graph.
generateGraph=False

# Choose options for the graph. Will only be used if generateGraph is set
#to True.

# Boolean flag to signal the choice of zooming in (or not) on the most
#important region
magnified=True

# Limits of the list of radii whose ground states will be computed. 
# Only used if the flag 'magnified' is set to True.
start=10
end=50

# Number of shots to be used in obtaining the expectation values using the
#CIRQ simulator
repetitions=1000

# Number of runs. If runs>1, the median over all the runs will be taken.
# Should be an odd number.
runs=1

if generateGraph:

    if (magnified):
        # Select desired window of the list of radii
        listOfRadii=listOfRadii[(start+1):end]
        listOfHamiltonians=listOfHamiltonians[start+1:end]
        
    # Get list of ground energies and states from exact diagonalization
    exactGroundEnergies,exactGroundStates=\
        groundStatesFromDiagonalization(listOfHamiltonians)
        
    # Get list of ground energies and states from optimizing the analytically
    #calculated energy
    optimizedGroundEnergies,optimizedGroundStates=\
        groundStatesFromTheoreticalOptimization(listOfHamiltonians,runs)
    
    # Get list of ground energies and states from optimizing the energy
    #obtained from simulating the circuit
    groundEnergiesVQE,groundStatesVQE=\
        groundStatesFromVQE(listOfHamiltonians,repetitions,runs)
            
    # Calculate overlaps between the ground state obtained in VQE and the one
    #obtained from exact diagonalization
    overlapList=[]
    for (i,groundStateVQE) in enumerate(groundStatesVQE):
        overlapList.append(calculateOverlap(groundStateVQE,exactGroundStates[i]))
    
    f,axs=plt.subplots(2,figsize=[8,8],constrained_layout=True,sharex=True)
    plt.suptitle('Bond Dissociation Curve of He-H$^+$ ('+str(repetitions)+\
                 ' repetitions in QEE)')
    
    axs[1].set_xlabel('Atomic Separation (pm)')
    axs[1].set_ylabel('State Overlap')
    
    axs[0].set_ylabel('Ground Energy (MJ mol$^{-1}$)')
        
    theoretical=axs[0].scatter(listOfRadii,optimizedGroundEnergies,marker='o',\
                               color='silver')
    
    experimental=axs[0].scatter(listOfRadii,groundEnergiesVQE,marker='x',\
                                color='b')
    
    exact=axs[0].plot(listOfRadii,exactGroundEnergies,color='r')
    
    axs[0].legend([theoretical,experimental,exact[0]],\
                  ['Theoretical QEE','Experimental QEE','Exact Diagonalization'])
        
    axs[1].plot(listOfRadii,overlapList)
    
    plt.show()