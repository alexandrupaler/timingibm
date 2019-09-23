# Experiments on the ibmq_16_melbourne chip

## Method
* Find the qubit with shortest T1
* Initialise qubit in |1>
* Execute a sequence of layers of CNOTs
* The last circuit operation is to measure the qubit for being in |1>

## Construction of circuits
* The SEQ layers include a single CX
* The PAR layers include 5 gates ['CX13_1', 'CX12_2', 'CX11_3', 'CX4_10', 'CX6_8'] 
By the topology of the chip the CNOTs in PAR layers are parallel.

## Assumption 

is that execution time for SEQ and PAR is approximately the same, 
thus at the end of the circuit, the probability of the selected T1 qubit being |1> is 
approximately equal between the scenarios.

In case the remote/cloud quantum computer would execute the PAR circuits sequentially 
then the execution times would have to be roughly 5 times longer, and the observed
probabilities for |1> would be very different. 
(e.g. when for PAR the prob is very close to zero, this would not be the case for SEQ) 


## Result

Shots, Layers, SEQ, PAR

8192,  4, 5062, 5044

8192,  6, 4440, 4238

8192,  8, 3636, 3629

8192, 10, 3148, 2958

8192, 12, 2588, 2721

8192, 14, 1961, 2112

8192, 21, 1336, 1445