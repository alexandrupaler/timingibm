from qiskit.result import Result
import matplotlib

import matplotlib.pyplot as plt
import numpy as np

from itertools import combinations

def get_bit_combinations(distance, ignore_qubits = []):
    nrqubits = 14
    tuples = list(combinations(range(nrqubits), distance))

    masks = []
    for i in range(len(ignore_qubits) + 1):
        masks += list(combinations(ignore_qubits, i))

    ret = []
    for t in tuples:

        if any(map(lambda x: x in ignore_qubits, t)):
            # if any qubit from the tuple should be ignored
            # then the entire tuple should be discarded
            # print("ignore", t)
            continue

        for m in masks:
            tt = list(t) + list(m)
            r = map(lambda x: 1 << x, tt)
            ret.append(r)

    return ret

def create_plot(results):
    # Data for plotting
    xaxis = list(map(hex, range(8192)))

    #get max value for non SIM
    maxv = -1
    for res in results.keys():
        if res != "SIM":
            for val in results[res]:
                maxv = max(maxv, val)

    print(maxv)

    plt.ylim(top = maxv * 1.2)

    for res in results.keys():
        if res != "SIM":
            plt.plot(xaxis, results[res], label=res)

    plt.legend()

    # fig, ax = plt.subplots()
    # ax.plot(xaxis, s)
    #
    # ax.set(xlabel='Basis', ylabel='Nr.',
    #        title='Histogram of measurements')
    # ax.grid()
    plt.savefig("test.png")

    # plt.show()

# def get_distance_basis(dist, nrqubits):
#     ret = []
#     for i in range(nrqubits):
#
#
#     return ret


def append_result(total_res, resname, dictionary):

    print("process {} with {} keys".format(resname, len(dictionary)))

    part_res = []

    wherein = 0
    total = 0
    for i in range(1<<14):
        val = 0
        if i in dictionary.keys():
            val = dictionary[i]
            wherein += 1
            total += val
        part_res.append(val)

    assert(total == 8192)
    assert(wherein == len(dictionary))

    total_res[resname] = part_res
    return total_res

def print_state(basis, mask_qubits):
    str = list("{0:014b}".format(basis))
    for x in mask_qubits:
        str[13 - x] = '*'

    print("    " + "".join(str))

def analysis_results(total_res, dists, mask_qubits=[]):

    for resname in total_res:
        if resname == "SIM":
            print(total_res[resname])
        else:
            print("\n" + resname)
            everythingshouldbe8192 = sum(total_res[resname])
            grandtotal = 0
            for dist in dists.keys():
                total = 0
                for d in dists[dist]:
                    x = total_res[resname][d]
                    # print_state(d, mask_qubits)
                    total += x

                print("DIST {} {}".format(dist, total))

                grandtotal += total
            print("GRAND {}/{}".format(grandtotal, everythingshouldbe8192))



def get_distance_keys(original, max_dist, ignore_qubits=[]):
    dist_key_basis = {}
    for dist in range(max_dist + 1):

        dist_key_basis[dist] = []

        bits_to_flip = get_bit_combinations(dist, ignore_qubits=ignore_qubits)

        for bf in bits_to_flip:
            val = original
            for b in bf:
                val = val ^ b
            dist_key_basis[dist].append(val)

    return dist_key_basis


def main():
    total_res = {}

    """
        ONE CONSTRUCT
    """
    # fname = "res_Tue Sep 24 14:55:10 2019.txt.experiment_two"
    # mask_qubits = [0, 7]

    """
        ONE DECAY
    """
    # fname = "res_Tue Sep 24 22:00:36 2019.txt"
    # mask_qubits = [0, 1]

    fname = "res_Mon Nov 25 00:50:09 2019.txt_1111_experiment"
    mask_qubits = [0, 7]

    """
        PLUS DECAY
    """
    # short ghz
    # fname = "res_Mon Nov 25 00:36:47 2019.txt_long_waiting_times"
    # mask_qubits = [0, 1, 6, 7, 8, 9, 10]
    # large ghz
    # fname = "res_Mon Nov 25 00:36:47 2019.txt_long_waiting_times"
    # mask_qubits = [0, 1]

    # 0b11111101111110
    # bits 0 and 7 are ZERO
    # all other bits are ONE
    # the perfect result will be read from the specified file
    perfect_result = 0x3F7E

    with open(fname, "r") as file:
        readnext = False
        resname = ""

        for line in file:

            if line.startswith(" ---"):
                resname = line.replace(" ---", "").strip()

            if line.startswith("{"):
                dict = eval(line)
                res = Result.from_dict(dict)

                for j in range(len(res.results)):
                    expcounts = {}
                    dict = res.results[j].data.counts.to_dict()
                    for key in dict.keys():
                        expcounts[int(key, 16)] = dict[key]

                    if resname == "SIM":
                        # this is the perfect basis vector from the simulator
                        perfect_result = list(expcounts.keys())[0]
                        print("the perfect result is", perfect_result, hex(perfect_result))
                        print_state(perfect_result, mask_qubits=mask_qubits)

                    circ_name = res.results[j].header.name
                    append_result(total_res, resname + circ_name, expcounts)

    dists = get_distance_keys(perfect_result, 3, ignore_qubits=mask_qubits)
    analysis_results(total_res, dists, mask_qubits=mask_qubits)

    # create_plot(total_res)


if __name__ == "__main__":
    main()