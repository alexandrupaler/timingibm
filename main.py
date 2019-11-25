import builtins as __builtin__

import time
import csv

import qiskit

from qiskit import IBMQ, execute
from qiskit.visualization import circuit_drawer

from qiskit.compiler import transpile
from qiskit.providers.jobstatus import JOB_FINAL_STATES

# coupling_map = [[1,0], [1,2], [2,3], [4,3], [4,10], [5,4], [5,6], [5,9], [6,8], [7,8], [9,8], [9,10], [11,10], [11,3], [11,12], [12,2], [13,1], [13,12]]

def print(*args, **kwargs):
    fconsole.write(''.join(map(lambda x : str(x) + ' ', args)))
    fconsole.write('\n')

    return __builtin__.print(*args, **kwargs)


def compute_cnot_time(filtered_cnots):
    # From: https://github.com/Qiskit/ibmq-device-information/tree/master/backends/melbourne/V1
    # A frame change (FC) is equivalent to applying a virtual Z-gate in software,
    # where Z(θ)=FC(-θ). Gaussian derivative (GD) and Gaussian flattop (GF) pulses
    # are defined with amplitude and angle parameters.
    #
    # There is an additional buffer of 20 ns after each GD or GF pulse.
    buffer = 20
    # All the GD have a gate time of 100 ns
    GD = 100
    # see https://github.com/Qiskit/ibmq-device-information/blob/master/backends/melbourne/V1/version_log.md
    # I copy pasted from there in the following file
    cnot, GF = select_value(max, "ibmq_16_times.txt",
                            name=0, val=1,
                            skip_first_line=False,
                            keys_from_list=filtered_cnots)

    print("->", cnot, GF)
    #
    time = GD + buffer + GF + buffer + GD + buffer + GF + buffer

    return cnot, time


def select_value(function_cmp, csv_file_name, name = 0, val = 1, skip_first_line = True, keys_from_list=None):
    values = [-1, 100000]

    value = [item for item in values if item not in [function_cmp(values)]][0]
    key = value

    with open(csv_file_name, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')

        #skip first row
        if skip_first_line:
            next(csvreader)

        # T1 is second column
        for row in csvreader:

            # skip the keys which are not allowed
            if keys_from_list is not None:
                if row[name] not in keys_from_list:
                    continue

            curr_value = float(row[val])
            value = function_cmp(value, curr_value)

            if value == curr_value:
                key = row[name]


    return key, value

def create_circuit_one(run_sequential_experiment = False, nr_layers = 1, list_of_cnots = [], t1_qubit = 0):
    qr = qiskit.QuantumRegister(14)
    cr = qiskit.ClassicalRegister(1)
    qc = qiskit.QuantumCircuit(qr, cr)

    # Flip the qubit with longest T1
    qc.x(qr[int(t1_qubit)])
    qc.barrier()

    barrier_u1_qubits = qr[0:7]

    if run_sequential_experiment:
        for i in range(nr_layers):
            qs = list_of_cnots[0].replace("CX", "").split("_")
            qc.cx(qr[int(qs[0])], qr[int(qs[1])])
            qc.barrier()
    else:
        for i in range(nr_layers):
            for gate in list_of_cnots:
                qs = gate.replace("CX", "").split("_")
                qc.cx(qr[int(qs[0])], qr[int(qs[1])])
                # qc.barrier()
            qc.barrier()

    qc.barrier()
    qc.measure(qr[int(t1_qubit)], cr[0])

    return qc

def create_circuit_one_create(run_sequential_experiment = False, list_of_cnots = []):
    qr = qiskit.QuantumRegister(14)
    cr = qiskit.ClassicalRegister(14)
    qc = qiskit.QuantumCircuit(qr, cr)

    ones = [4, 5, 6, 11, 12, 13]
    for q in ones:
        qc.x(qr[q])
    qc.barrier()

    for gate in list_of_cnots:
        qs = gate.replace("CX", "").split("_")
        qc.cx(qr[int(qs[0])], qr[int(qs[1])])

        if run_sequential_experiment:
            qc.barrier()

    qc.barrier()
    # for i in [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13]:
    for i in range(14):
        qc.measure(qr[i], cr[i])

    return qc


def create_circuit_one_decay(nr_layers = 1, circuit_name = None):
    qr = qiskit.QuantumRegister(14)
    cr = qiskit.ClassicalRegister(14)
    qc = qiskit.QuantumCircuit(qr, cr, name = circuit_name)

    clocks = [0, 1]
    for q in clocks:
        qc.x(qr[q])
    qc.barrier()

    ones = []
    for x in range(14):
        if x not in clocks:
            ones.append(x)

    # Create the 11111101111110 state
    for q in ones:
        qc.x(qr[q])
    qc.barrier()

    for i in range(nr_layers):
        # Let the qubits do nothing
        qc.cx(qr[1], qr[0])
        qc.barrier()

    qc.barrier()
    for i in range(2, 14):
        qc.measure(qr[i], cr[i])

    return qc

def create_circuit_plus_decay(nr_layers = 1, circuit_name = None):
    qr = qiskit.QuantumRegister(14)
    cr = qiskit.ClassicalRegister(14)
    qc = qiskit.QuantumCircuit(qr, cr, name = circuit_name)

    clocks = [0, 1]
    for q in clocks:
        qc.x(qr[q])
    qc.barrier()

    ones = []
    for x in range(14):
        if x not in clocks:
            ones.append(x)

    # Construct the logical_+ state
    # qc.h(qr[ones[0]])
    # for q in range(1, len(ones)):
    #     qc.cx(qr[ones[0]], qr[ones[q]])
    # qc.barrier()

    qc.h(qr[13])
    qc.cx(qr[13], qr[12])
    qc.cx(qr[12], qr[2])
    qc.cx(qr[2], qr[3])

    qc.h(qr[3])
    qc.h(qr[4])
    qc.h(qr[11])

    qc.cx(qr[4], qr[3])
    qc.cx(qr[11], qr[3])

    qc.h(qr[11])

    qc.h(qr[5])
    qc.cx(qr[5], qr[4])
    qc.h(qr[5])

    qc.h(qr[4])
    qc.h(qr[3])
    # |13 12 2 3 4 11 5>
    # mask = [0, 1] + unused
    # mask = [0, 1, 6, 7, 8, 9, 10]


    # qc.cx(qr[11], qr[10])
    # qc.cx(qr[5], qr[9])
    # qc.cx(qr[5], qr[6])
    # qc.cx(qr[6], qr[8])
    #
    # qc.h(qr[8])
    # qc.h(qr[7])
    # qc.cx(qr[7], qr[8])
    # qc.h(qr[8])
    # qc.h(qr[7])

    qc.barrier()


    for i in range(nr_layers):
        # Let the qubits do nothing
        qc.cx(qr[1], qr[0])
        qc.barrier()

    qc.barrier()

    # # Destruct the logical_+ state
    # for q in range(1, len(ones)):
    #     qc.cx(qr[ones[0]], qr[ones[q]])
    # qc.h(qr[ones[0]])
    #
    # qc.barrier()

    # qc.h(qr[8])
    # qc.h(qr[7])
    # qc.cx(qr[7], qr[8])
    # qc.h(qr[8])
    # qc.h(qr[7])
    #
    # qc.cx(qr[6], qr[8])
    # qc.cx(qr[5], qr[6])
    # qc.cx(qr[5], qr[9])
    # qc.cx(qr[11], qr[10])
    #
    qc.h(qr[5])

    qc.h(qr[4])
    qc.h(qr[3])
    qc.cx(qr[5], qr[4])

    qc.h(qr[11])

    qc.h(qr[5])

    qc.cx(qr[11], qr[3])
    qc.cx(qr[4], qr[3])

    qc.h(qr[3])
    qc.h(qr[4])
    qc.h(qr[11])

    qc.cx(qr[2], qr[3])
    qc.cx(qr[12], qr[2])
    qc.cx(qr[13], qr[12])

    qc.h(qr[13])

    qc.barrier()


    for i in range(2, 14):
        qc.measure(qr[i], cr[i])

    return qc

def run_on_ibm(q_circuit, backend):
    """
        I am setting optimization_level=0 in the hope that gates are not decomposed or optimised between any pair of barriers
    """
    qobj = qiskit.transpile(q_circuit, backend, optimization_level=0)
    #
    job = execute(qobj, backend=backend, shots=8192)

    start_time = time.time()
    job_status = job.status()
    while job_status not in JOB_FINAL_STATES:
        print(f'Status @ {time.time()-start_time:0.0f} s: {job_status.name},'
              f' est. queue position: {job.queue_position()}')
        time.sleep(10)
        job_status = job.status()

    return job.result()

def experiment_one_state_create(back_ibm, back_sim):
    """
    The goal of this experiment is to investigate if the same set of CNOTs introduces
    more errors when executed in parallel or sequentially.

    Therefore, two logical qubits are considered
    q_1 = |111000>
    q_2 = |000111>

    The states have been chosen based on the topology of ibmq_16_melbourne with 6 qubits on the top row,
    and 6 qubits on the lower row. Qubits Q0 (top-left) and Q7 (bottom-right) are not used for this experiment.


    Six CNOTs are executed between q_1 and q_2. The individual CNOTs are executed in SEQ and PAR.
    SEQ has barriers between the CNOTs, while PAR assumes that all six CNOTs parallel.

    Note: It could be that X errors introduced by the CNOTs cancel out or are cascading.
    Nevertheless, it is interesting to observe what the codeword at the end of the experiment is.

    At the end the qubits are measured. The perfect state should be
    q_1=|111111>
    q_2=|111111>

    The scripts in resanalysis.py can analyse the measurement results and compute
    * the incidence of measuring the perfect |q_1q_2> = |1111111_111111_> where _ (underscore) is Q0 and Q7
    * the cumulated incidence of distances D (e.g. 1) basis from the perfect |q_1q_2>
    """

    parallel_cnots = ['CX13_1', 'CX12_2', 'CX11_3', 'CX4_10', 'CX5_9', 'CX6_8']

    qc_seq = create_circuit_one_create(run_sequential_experiment=True,
                                       list_of_cnots=parallel_cnots)

    #
    res_seqs = run_on_ibm(qc_seq, back_sim)
    print(' --- SIM')
    print(res_seqs.to_dict())
    #
    res_seq = run_on_ibm(qc_seq, back_ibm)
    print(' --- SEQ')
    print(res_seq.to_dict())

    qc_par = create_circuit_one_create(run_sequential_experiment=False,
                                       list_of_cnots=parallel_cnots)
    res_par = run_on_ibm(qc_par, back_ibm)
    print(' --- PAR')
    print(res_par.to_dict())

    # Print the ASCII circuits
    print(circuit_drawer(qc_seq, line_length=-1))
    print(circuit_drawer(qc_par, line_length=-1))

def experiment_one_state_decay(back_ibm, back_sim):
    """
        Create a circuit with variable number of layers that do nothing
        In order to not perform any operation, but also to speed up the experiment
        two qubits are chosen to perform CNOT sequentially between barriers
        Each circuit layer has the duration of the CNOT between those qubits

        The question asked by this sequence of experiments is:
            "How long does it take until the logical |1> state on all the other qubits
            becomes a logical |0> state?"

        There are 12 qubits in the codeword, thus the question is:

            "How does the probability of distance less than 5 errors depend with the number
            CNOT layers in the circuit?"

            "In terms of CNOT times: how long does it take until the codeword is affected beyond recognition by
            distance 6 or higher errors?"
    """

    qc_seq = create_circuit_one_decay(nr_layers=1)
    # res_seqs = run_on_ibm(qc_seq, back_sim)
    # print(' --- SIM')
    # print(res_seqs.to_dict())
    #
    # for i in range(1, 11):
    #     qc_seq = create_circuit_one_decay(nr_layers=i * 2)
    #     res_seqs = run_on_ibm(qc_seq, back_ibm)
    #     print(' --- LAY{}'.format(2 * i))
    #     print(res_seqs.to_dict())
    #
    # print(circuit_drawer(qc_seq, line_length=-1))

    res_seqs = run_on_ibm(qc_seq, back_sim)
    print(' --- SIM')
    print(res_seqs.to_dict())

    print(' --- HRD')

    qc_seqs = []
    for i in range(1, 10):
        nr_layers = 2 * i
        cname = " --- LAY_{}".format(nr_layers)
        qc_seq = create_circuit_one_decay(nr_layers=nr_layers,
                                           circuit_name=cname)
        qc_seqs.append(qc_seq)
        # res_seqs = run_on_ibm(qc_seq, back_ibm)

    res_seqs = run_on_ibm(qc_seqs, back_ibm)
    print(res_seqs.to_dict())

    print(circuit_drawer(qc_seq, line_length=-1))



def experiment_plus_state_decay(back_ibm, back_sim):
        """
            Create a circuit with variable number of layers that do nothing
            In order to not perform any operation, but also to speed up the experiment
            two qubits are chosen to perform CNOT sequentially between barriers
            Each circuit layer has the duration of the CNOT between those qubits

            The question asked by this sequence of experiments is:
                "How long does it take until the logical |1> state on all the other qubits
                becomes a logical |0> state?"

            There are 12 qubits in the codeword, thus the question is:

                "How does the probability of distance less than 5 errors depend with the number
                CNOT layers in the circuit?"

                "In terms of CNOT times: how long does it take until the codeword is affected beyond recognition by
                distance 6 or higher errors?"
        """

        qc_seq = create_circuit_plus_decay(nr_layers = 1,
                                           circuit_name = "simulation")
        res_seqs = run_on_ibm(qc_seq, back_sim)
        print(' --- SIM')
        print(res_seqs.to_dict())

        print(' --- HRD')

        qc_seqs = []
        for i in range(1, 10):
            # short with long times
            # nr_layers = 30 * i
            # long with short times
            nr_layers = 2 * i
            cname = " --- LAY_{}".format(nr_layers)

            qc_seq = create_circuit_plus_decay(nr_layers = nr_layers,
                                               circuit_name = cname)
            qc_seqs.append(qc_seq)
            # res_seqs = run_on_ibm(qc_seq, back_ibm)

        res_seqs = run_on_ibm(qc_seqs, back_ibm)
        print(res_seqs.to_dict())

        print(circuit_drawer(qc_seq, line_length=-1))

def main():

    # IBMQ.save_account('9f4a49cb841ae1e3ce8b1de1d6a7860e1c153019d274d3fb1c270e9da6f17e420c3a37c4ceb737d90d7c80908f642f9c101aabba9601a7753d4b7f30a156d467')

    IBMQ.enable_account(
            "9f4a49cb841ae1e3ce8b1de1d6a7860e1c153019d274d3fb1c270e9da6f17e420c3a37c4ceb737d90d7c80908f642f9c101aabba9601a7753d4b7f30a156d467")

    back_ibm = IBMQ.get_provider(hub='ibm-q', group='open', project='main').backends(name='ibmq_16_melbourne')[0]
    back_sim = IBMQ.get_provider(hub='ibm-q', group='open', project='main').backends(name='ibmq_qasm_simulator')[0]

    experiment_plus_state_decay(back_ibm, back_sim)

    # experiment_one_state_decay(back_ibm, back_sim)

    # qubit, t1_time = select_value(min, "ibmq_16_melbourne.csv", name=0, val=1)
    # # microseconds to nanoseconds
    # t1_time *= 1000
    # print("qubit, t1_time", qubit, t1_time)
    #
    # """
    # Create the list of CNOTs from a layer
    # """
    # qubit = qubit.replace("Q", "")
    # # number of max parallel gates
    # # I am choosing the vertical cnots between the horizontal lines
    # parallel_cnots = ['CX13_1', 'CX12_2', 'CX11_3', 'CX4_10', 'CX5_9', 'CX6_8']
    # filtered_cnots = []
    # # Do not apply any gates on the T1 qubit
    # for gate in parallel_cnots:
    #     qubits = gate.replace("CX", "").split("_")
    #     if qubit not in qubits:
    #         filtered_cnots.append(gate)
    #
    #
    # """
    # Maximum duration of a CNOT
    # """
    # selected_cnot, max_cx_time = compute_cnot_time(filtered_cnots)
    # print("cnot time", max_cx_time)
    #
    # cnots_in_a_layer = len(filtered_cnots)
    # print(filtered_cnots)
    #
    # # when things go strange, multiply with a factor
    # factor = 0.4
    #
    # # how many sequential cx gates are min to achieve t1_time?
    # # assuming that all gates in a layer have the longest cx_time
    # # Thus, after ratio layers the state of qubit for which T1 is considered should be |0>
    # nr_seq_cx = round((t1_time * factor)/ max_cx_time)
    #
    # #     FORCE TO 1
    # # nr_seq_cx = 1
    #
    # print(nr_seq_cx)
    # # qc_seq = create_circuit_one(run_sequential_experiment=True,
    # #                         nr_layers = nr_seq_cx,
    # #                         list_of_cnots = [selected_cnot],
    # #                         t1_qubit=qubit)
    # qc_seq = create_circuit_two(run_sequential_experiment=True,
    #                         list_of_cnots = parallel_cnots)
    # #
    # res_seqs = run_on_ibm(qc_seq, back_sim)
    # print(' --- SIM')
    # print(res_seqs.to_dict())
    # #
    # res_seq = run_on_ibm(qc_seq, back_ibm)
    # print(' --- SEQ')
    # print(res_seq.to_dict())
    #
    # # I am placing the sequential CX in layers of gates which I assume are executed in parallel
    # # If the gates within the layer are executed in parallel then
    # # the state of the T1 qubit should be prob < 1 in a state diff from |0>
    # # otherwise with prob almost one the qubit is in state |0>
    # # qc_par = create_circuit_one(run_sequential_experiment=False,
    # #                             nr_layers = nr_seq_cx,
    # #                             list_of_cnots = filtered_cnots,
    # #                             t1_qubit=qubit)
    #
    # qc_par = create_circuit_two(run_sequential_experiment=False,
    #                             list_of_cnots=parallel_cnots)
    # res_par = run_on_ibm(qc_par, back_ibm)
    # print(' --- PAR')
    # print(res_par.to_dict())
    #
    # # Print the ASCII circuits
    # print(circuit_drawer(qc_seq, line_length=-1))
    # print(circuit_drawer(qc_par, line_length=-1))

    return

if __name__ == "__main__":
    with open("res_{}.txt".format(time.asctime()), "w") as fconsole:
        # read the original script before doing anything else
        script_lines = []
        with open("main.py") as src:
            script_lines = src.readlines()

        # run main loop
        main()

        # write the script to the res
        for i in range(10):
            fconsole.write('\n')
        fconsole.write('# the code for this experiment is below\n')
        fconsole.writelines(script_lines)