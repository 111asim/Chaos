
# import numpy as np
# from qiskit import QuantumCircuit, transpile,Aer,execute
# from qiskit.quantum_info import Kraus, SuperOp
from qiskit_aer import AerSimulator,qasm_simulator
from qiskit.tools.visualization import plot_histogram
from qiskit.quantum_info import DensityMatrix, partial_trace,random_statevector
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from qiskit.visualization import plot_state_city, plot_state_qsphere
import itertools

def get_cyclic_permutation(permutation):
    '''
    input: arbitary permutation of tokens in complete graph:list
    output: cyclic permutation equivalent of the input:list of lists
    '''
    pi=permutation
    cyclic_permutation=[]
    n=len(pi)
    not_in_cycle=list(range(n))
    while len(not_in_cycle)!=0:
        cycle=[]
        start=not_in_cycle[0]
        current=start
        cycle.append(start)
        not_in_cycle.remove(start)
        final=-1
        while final != start:
            current=permutation[current]
            cycle.append(current)
            if(current in not_in_cycle):
                not_in_cycle.remove(current)
            final=current
        cyclic_permutation.append(cycle)
    cyclic=[]
    for i in cyclic_permutation:
        length=len(i)
        if(length>2):
            cyclic.append(i[:-1])
    return(cyclic)
def Get_Routing_via_matching(cyclic_permutation):
    '''
    input:cyclic permutations: list of lists
    output: 2 step swaps that is equivalent to routing via matching 2 list of list(list of swaps)
    '''
    layer_1=[]
    layer_2=[]
    for cycle in  cyclic_permutation:
        length=len(cycle)
        for j in range(int(len(cycle)/2)):
            layer_1.append([cycle[j],cycle[length-j-1]])
            if(cycle[j+1]!=cycle[length-j-1]):
                layer_2.append([cycle[j+1],cycle[length-j-1]])
    return(layer_1,layer_2)

def old_causaly_cover_nodes(M,i,j):##edges=M
    '''
    input: M=random matching for the complete graph
           i=initial node n>=i,j>=0 int
           j=final node 
    output:boolean:true: iand j are causaly cover
                   false: i and j are not causaly covered
    '''
    path=[i]
    cover=False
    for matching in M: ##if you think about it M is nothing but random matching in our complete graph
        for i in matching:
            current=path[-1]
            if(i[0]==current):
                path.append(i[1])
                if(i[1]==j):
                    cover=True
                break
            elif(i[1]==current):
                path.append(i[0])
                if(i[0]==j):
                    cover=True
                break
        if(cover==True):
            # print("path for causaly cover=",path)
            return cover
    return cover
def causaly_cover_nodes(M,i,j):##edges=M
    '''
    input: M=random matching for the complete graph
           i=initial node n>=i,j>=0 int
           j=final node 
    output:boolean:true: iand j are causaly cover
                   false: i and j are not causaly covered
    '''
    covered_set=[i]
    cover=False
    for matching in M: ##if you think about it M is nothing but random matching in our complete graph
        for k in matching:
            new_covered=[]
            for l in covered_set:
                current=l
                if(k[0]==current):
                    new_covered.append(k[1])
                    if(k[1]==j):
                        cover=True
                        break
                elif(k[1]==current):
                    new_covered.append(k[0])
                    if(k[0]==j):
                        cover=True
                        break
            covered_set=covered_set+new_covered
            if(cover==True):
                # print("path for causaly cover=",path)
                return cover
    return cover

def causaly_cover_graph(M,n):
    nodes=list(range(n))
    ##get all the possible pair combination for nodes:
    permutation=list(itertools.permutations(nodes,2))
    for pairs in permutation:
        is_causaly_cover=causaly_cover_nodes(M,pairs[0],pairs[1])
        if(is_causaly_cover== False):
            return(False,pairs)
    return(True,(-1,-1))
def create_M(pi):
    ''' 
    input: random/ arbitary permutations
    output:return all the matching generated for this permutations using 2 step routing algorithm
    '''
    ##pass all these permutation to get two step routing 
    M=[]
    for i in range(len(pi)):
        cyclic_permutation=get_cyclic_permutation(pi[i])
        O,E=Get_Routing_via_matching(cyclic_permutation)
        M.append(O)
        M.append(E)
    return M


def initialize_zero(qc,n):
    for i in range(n):
        qc.initialize([1,0],i)
    return qc

def initialize_R_product_state(q,n):
    for i in range(n):
        q.initialize(random_statevector(2).data,i)
    return(q)

def random_clifford(qc,n):
    qubits_list=list(range(n))
    for i in range(500):
        choice=random.choices(['s','h','cnot','i'])[0]
        qbit=random.choices(qubits_list)[0]
        if choice=='s':
            qc.s(qbit)
        elif choice=='h':
            qc.h(qbit)
        elif choice=='i':
            True
        elif choice=='cnot':
            c=qubits_list.copy()
            c.remove(qbit)
            control_bit=random.choices(c)[0]
            order=[control_bit,qbit]
            random.shuffle([control_bit,qbit]) ## randomize the cont order
            qc.cx(order[0],order[1])
    return qc

def random_clifford_with_depth(qc,n,depth):
    qubits_list=list(range(n))
    gates=[]
    for i in range(depth):
        un_used_qubits=qubits_list.copy()
        gates_added={}
        while len(un_used_qubits)>=1:
            qbit=random.choice(un_used_qubits)
            choice=random.choices(['s','h','cnot'])[0]
            un_used_qubits.remove(qbit)
            if choice=='s':
                qc.s(qbit)
                gates_added[qbit]='s'
            elif choice=='h':
                qc.h(qbit)
                gates_added[qbit]='h'
            elif choice=='cnot':
                if(un_used_qubits!=[]):
                    control_bit=random.choices(un_used_qubits)[0]
                    order=[control_bit,qbit]
                    random.shuffle(order)
                    qc.cx(order[0],order[1])
                    gates_added[str(order)]='cnot'

        gates.append(gates_added)
    return qc,gates

def random_Cnot_permutation(qc,n):
    pi=[]
    nodes=list(range(n))
    for i in range(int(math.log(n,2))):
        pi.append(list(random.sample(nodes,n)))
    M=create_M(pi)
    ## now add the cnot gate layers according to the matching. 
    for i in M:
        for cnots in i:
            qc.cx(cnots[0],cnots[1])
        ## after each layer of cnots add random S,H gates to all other qubits
        for i in range(n):
            choice=random.choices(['s','h','i'])[0]
            if choice=='s':
                qc.s(i)
            elif choice=='h':
                qc.h(i)
            elif choice=='i':
                True
    return qc

def t_layer(qc,n):
    for i in range(n):
        qc.t(i)
    return qc

def get_r_tildae_qiskit(qc,partition):
    rho= DensityMatrix.from_instruction(qc) ##get the state vector of the system
    # print(state)
    rho_a=partial_trace(state=rho,qargs=partition)  ##rho_a is the reduced density matrix for the first a subsystem A(as the paper has divided the whole system into A and B subsystem)
    eigenvalues=np.linalg.eigvals(rho_a.data) ##entanglement spectrum is the set of the eigen values of the reduced density matrix
    ES=sorted(eigenvalues,reverse=True)
    ES_real=[i.real for i in ES ]
    lambd=ES_real
    r_k_tilda=[]
    for i in range(1,len(lambd)-1):
        del_k=lambd[i-1]-lambd[i]
        del_k_P1=lambd[i]-lambd[i+1]
        r_k_tilda.append(min(del_k,del_k_P1)/max(del_k,del_k_P1))
    r_mean=np.mean(r_k_tilda)
    return r_mean

def rdm_qiskit(qc,partition):
    rho= DensityMatrix.from_instruction(qc) ##get the state vector of the system
    # print(state)
    rho_a=partial_trace(state=rho,qargs=partition)  ##rho_a is the reduced density matrix for the first a subsystem A(as the paper has divided the whole system into A and B subsystem)
    return rho_a

def brick_1D_add_4_layer(qc,n): 
    layer=0
    for i in range(0,n,2):
        pair=[i,i+1]
        random.shuffle(pair)
        qc.cx(pair[0],pair[1])
    for i in range(n):
        choice=random.choices(['s','h','i'])[0]
        if choice=='s':
            qc.s(i)
        elif choice=='h':
            qc.h(i)
        elif choice=='i':
            True
    for i in range(1,n-1,2):
        pair=[i,i+1]
        random.shuffle(pair)
        qc.cx(pair[0],pair[1])
    # qc.cnot(0,n-1)
    for i in range(n):
        choice=random.choices(['s','h','i'])[0]
        if choice=='s':
            qc.s(i)
        elif choice=='h':
            qc.h(i)
        elif choice=='i':
            True
    return qc

def rdm_blue_qubit(qc,partition=0):
    n=qc.num_qubits
    if (partition==0):
        partition=list(range(int(n/2)))
    ## this loop for arranging the qubits incase the qubits for RDM is not sequencial
    qc_copy=qc.copy()
    partition.sort()
    for i in range(len(partition)):
        if(i!=partition[i]):
            qc_copy.swap(i,partition[i])
    qc_copy.measure_all()
    ##simulate using blue qubit
    bq=bluequbit.init("bnK9f3RMNBjj2Fi3BqyBbcyFLHS8pLFc")
    result=bq.run(qc_copy,job_name="testing_1")
    statevector=result.get_statevector()
    state_matrix=statevector.reshape(2**int(n/2),2**int(n/2))
    rdm=np.matmul(state_matrix,np.conj(state_matrix).T)
    return rdm

def get_r_tildae_rdm(rdm):
    eigenvalues=np.linalg.eigvals(rdm) ##entanglement spectrum is the set of the eigen values of the reduced density matrix
    ES=sorted(eigenvalues,reverse=True)
    ES_real=[i.real for i in ES ]
    lambd=ES_real
    r_k_tilda=[]
    for i in range(1,len(lambd)-1):
        del_k=lambd[i-1]-lambd[i]
        del_k_P1=lambd[i]-lambd[i+1]
        r_k_tilda.append(min(del_k,del_k_P1)/max(del_k,del_k_P1))
    r_mean=np.mean(r_k_tilda)
    return r_mean

def get_r_tildae(qc,partition,simulator):
    n=qc.num_qubits
    if (partition==0):
        partition=list(range(int(n/2)))
    if(simulator=='qiskit'):
        rdm=rdm_qiskit(qc,partition)
        r_tildae=get_r_tildae_rdm(rdm)
    elif(simulator=='bq'):
        rdm=rdm_blue_qubit(qc,partition)
        r_tildae=get_r_tildae_rdm(rdm)
    return(r_tildae,rdm)