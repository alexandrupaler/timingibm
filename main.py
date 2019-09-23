import time
import csv

import qiskit

from qiskit import IBMQ, execute
from qiskit.compiler import transpile
from qiskit.providers.jobstatus import JOB_FINAL_STATES

coupling_map = [[1,0], [1,2], [2,3], [4,3], [4,10], [5,4], [5,6], [5,9], [6,8], [7,8], [9,8], [9,10], [11,10], [11,3], [11,12], [12,2], [13,1], [13,12]]

def get_maximum_gf():

    value = -1
    with open("ibmq_16_times.txt", newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')

        next(csvreader)



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


def create_circuit(run_sequential_experiment = False, nr_layers = 1, list_of_cnots = [], t1_qubit = 0):
    qr = qiskit.QuantumRegister(14)
    cr = qiskit.ClassicalRegister(1)
    qc = qiskit.QuantumCircuit(qr, cr)

    # Flip the qubit with longest T1
    qc.x(qr[int(t1_qubit)])
    qc.barrier()

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


def run_on_ibm(q_circuit):
    backend = IBMQ.get_provider(hub='ibm-q', group='open', project='main').backends(name='ibmq_16_melbourne')[0]
    # backend = IBMQ.get_provider(hub='ibm-q', group='open', project='main').backends(name='ibmq_qasm_simulator')[0]
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

#
# # IBMQ.save_account('9f4a49cb841ae1e3ce8b1de1d6a7860e1c153019d274d3fb1c270e9da6f17e420c3a37c4ceb737d90d7c80908f642f9c101aabba9601a7753d4b7f30a156d467')
#
IBMQ.enable_account(
        "9f4a49cb841ae1e3ce8b1de1d6a7860e1c153019d274d3fb1c270e9da6f17e420c3a37c4ceb737d90d7c80908f642f9c101aabba9601a7753d4b7f30a156d467")

qubit, t1_time = select_value(min, "ibmq_16_melbourne.csv", name=0, val=1)
# microseconds to nanoseconds
t1_time *= 1000
print("qubit, t1_time", qubit, t1_time)

"""
Create the list of CNOTs from a layer
"""
qubit = qubit.replace("Q", "")
# number of max parallel gates
# I am choosing the vertical cnots between the horizontal lines
parallel_cnots = ['CX13_1', 'CX12_2', 'CX11_3', 'CX4_10', 'CX5_9', 'CX6_8']
filtered_cnots = []
# Do not apply any gates on the T1 qubit
for gate in parallel_cnots:
    qubits = gate.replace("CX", "").split("_")
    if qubit not in qubits:
        filtered_cnots.append(gate)


"""
Maximum duration of a CNOT
"""
selected_cnot, max_cx_time = compute_cnot_time(filtered_cnots)
print("cnot time", max_cx_time)

cnots_in_a_layer = len(filtered_cnots)
print(filtered_cnots)

# when things go strange, multiply with a factor
factor = 2.0

# how many sequential cx gates are min to achieve t1_time?
# assuming that all gates in a layer have the longest cx_time
# Thus, after ratio layers the state of qubit for which T1 is considered should be |0>
nr_seq_cx = round((t1_time * factor)/ max_cx_time)
print(nr_seq_cx)
qc_seq = create_circuit(run_sequential_experiment=True,
                        nr_layers = nr_seq_cx,
                        list_of_cnots = [selected_cnot],
                        t1_qubit=qubit)
# print(qc_seq)
res_seq = run_on_ibm(qc_seq)
print(' --- SEQ')
print(res_seq)

# I am placing the sequential CX in layers of gates which I assume are executed in parallel
# If the gates within the layer are executed in parallel then
# the state of the T1 qubit should be prob < 1 in a state diff from |0>
# otherwise with prob almost one the qubit is in state |0>
max_nr_layers = round(nr_seq_cx / cnots_in_a_layer)
# qc_par = create_circuit(run_sequential_experiment=False, nr_layers = max_nr_layers, list_of_cnots = filtered_cnots)
qc_par = create_circuit(run_sequential_experiment=False,
                        nr_layers = nr_seq_cx,
                        list_of_cnots = filtered_cnots,
                        t1_qubit=qubit)
# print(qc_par)
res_par = run_on_ibm(qc_par)
print(' --- PAR')
print(res_par)

with open("res_{}.txt".format(time.asctime()), "w") as f:
    f.write("qubit {}, t1_time {}\n".format(qubit, t1_time))
    f.write("cnot time {}\n".format(max_cx_time))
    f.write('\n')
    f.write(str(res_seq))
    f.write('\n')
    f.write(str(res_par))
    for i in range(10):
        f.write('\n')
    f.write('# the code for this experiment is below\n')
    with open("main.py") as src:
        f.writelines(src.readlines())
