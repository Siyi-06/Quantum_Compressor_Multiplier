import matplotlib.pyplot as plt
import math
import numpy as np
from Startegy_Dynamic_Programming import dp
from fontTools.cu2qu.cu2qu import MAX_N


def Our_Design(N, max_m):
    # Parameters
    # N = 100  # Size of the multiplier (N x N)
    # max_m = 15  # Maximum value of m for m:k compressors
    max_layers = N  # Maximum number of layers in the Wallace tree, not include the last Carry Propagation Adder
    max_columns = 2 * N - 1  # Maximum columns in Wallace Tree
    P_value = [[0 for l in range(max_layers)] for j in range(max_columns)]
    # print(P_value)

    # First Layer: Initialize the partial product count
    for j in range(N):
        P_value[j][0] = j + 1
        # print("!", j, j + 1)
        if j != N - 1:
            P_value[max_columns - 1 - j][0] = j + 1
            # print("!",max_columns -1- j, j + 1)
    # print(P_value)

    # Cost of First Layer
    TC = 4 * N * N
    TD = N
    QC = 2 * N + N * N
    # print("Cost of First Layer:")
    # print("TC=",TC,"TD=",TD,"QC=",QC)

    # Compressor constraints and partial product updates
    flag_value = 3
    add_size = 100
    m_largest = 100
    m = 10000
    k = 10000
    for l in range(max_layers):
        # flag for check last layer
        m_largest = 0
        flag_value = 0
        add_size = 0
        for j in range(max_columns):
            if add_size == 0:
                if P_value[j][l] == 2:
                    add_size = max_columns - j
                    break
        for j in range(max_columns):
            if P_value[j][l] > 2:
                flag_value = 1
        # print("flag!!!",flag_value)
        if flag_value == 0:
            # print("max_l=",l)#maximum level 3 layers, not include CPA
            # print("add_size=", add_size)
            break
        for j in range(max_columns):
            # Apply compressors and update the partial products
            if P_value[j][l] > 1:
                m = P_value[j][l]  # use the largest compressor always
                k = math.floor(math.log2(m)) + 1
            if P_value[j][l] < 2:
                m = 0
                k = 0
                # print("P_value[", j, ",", l, "]", P_value[j][l], "m", m, "k", k)
                # print("No compressor")
                continue
            if l == 0:
                if P_value[j][l] > max_m:
                    m = max_m
                    k = math.floor(math.log2(m)) + 1
                    # print("YES")
                # print("P_value[", j, ",", l, "]", P_value[j][l], "m", m, "k", k)
            if l > 0:
                if P_value[j][l] > 2:
                    m = 3
                    k = 2
                if P_value[j][l] == 2:
                    m = 2
                    k = 2
                if P_value[j][l] < 2:
                    # print("No Compressor Here!")
                    m = 0
                    k = 0
                    continue
                # print("P_value[",j,",",l,"]",P_value[j][l],"m",m,"k",k)
                ###Update Cost
            # print("m_largest", m_largest)
            if m_largest < m:
                m_largest = m
                # print("m_largest",m_largest)
            if m > k:
                TC = TC + 4 * (m - k)
                QC = QC + (m - k)
            if m == k:
                TC = TC + 4
                QC = QC + 1
                # print("!!!B2")

                # Update the partial product count for column j after applying the compressor
                # Reduce the partial products at column j
            if l < (max_layers - 1):
                # Reduce the partial products at column j
                P_value[j][l + 1] = P_value[j][l] - m + 1 + P_value[j][l + 1]
                # print("Compression")
                # print("P_value[", j, ",", l, "]", P_value[j][l])
                # print("P_value[", j, ",", l+1, "]", P_value[j][l + 1])

                # Propagate carry-out bits to the next columns (j+1 to j+k-1)
                for ii in range(1, k):
                    if j + ii < 2 * N - 1:
                        # print("Before carry propagation")
                        # print("P_value[", j + ii, ",", l + 1, "]", P_value[j + ii][l + 1])
                        P_value[j + ii][l + 1] = P_value[j + ii][l + 1] + 1
                        # print("carry propagation")
                        # print("P_value[", j+ii, ",", l+1 , "]", P_value[j+ii][l+1])
        # Update Cost
        # print("Previous TD", TD)
        TD = TD + math.ceil(math.log(m_largest, 3)) + 2 * math.floor(math.log(m_largest, 2)) - 2
        # print("m_largest", m_largest)
        # print("TD",TD)
    # print(P_value)
    # print("Cost", "TC=",TC,"TD=",TD,"QC=",QC)
    # TC,TD,QC,add_size,max_l
    return TC, TD, QC, add_size, l


def BestRate_M(N):
    k = math.floor(math.log2(N)) + 1
    m = math.floor(math.pow(2, k) - 1)
    if m > N:
        k = k - 1
        m = math.floor(math.pow(2, k) - 1)
    if m == 2:
        k = 2
    if N == 2:
        m = 2
        k = 2
    if m < 2:
        m = 0
        k = 0
    # print("Best Rate M", m)
    return m


def Best_M(N, metric):
    flag = [math.inf, math.inf, math.inf, math.inf, math.inf]
    Best_M = 0
    for M in range(3, N + 1):
        # print(N,M)
        value = Our_Design(N, M)
        if value[metric] < flag[metric]:
            flag[metric] = value[metric]
            Best_M = M
    print('Best_M', Best_M)
    return Best_M


def cuccaro_qc(result):
    TC, TD, QC, add_size, l = result
    n = add_size
    TC += (2 * n - 1) * 7
    TD += (2 * n - 1) * 3
    QC += 2
    print(TC, TD, QC)
    return TC, TD, QC, add_size, l


def optimal_td(result):
    TC, TD, QC, add_size, l = result
    n = add_size
    TC += 2 * n * math.log2(n)
    TD += math.log2(n) + 1
    QC += n + n * math.log2(n) + math.ceil(math.log2(n)) + 2 - 2 * n
    print(TC, TD, QC)
    return TC, TD, QC, add_size, l


def gidney_tc(result):
    TC, TD, QC, add_size, l = result
    n = add_size
    TC += (n - 1) * 4
    TD += n
    QC += n - 1
    print(TC, TD, QC)
    return TC, TD, QC, add_size, l


if __name__ == "__main__":
    max_N = 39
    add_size_fns = {'Optimal TD': optimal_td,
                    'Optimal TC': gidney_tc}

    for METRIC_INDEX in [0]:
        x = []
        dp_td_y = []
        Largest_M_Mul_y = [[] for _ in range(len(add_size_fns))]
        Best_Rate_M_Mul_y = [[] for _ in range(len(add_size_fns))]
        Best_M_Mul_y = [[] for _ in range(len(add_size_fns))]
        Previous_Work_y = [[] for _ in range(len(add_size_fns))]

        METRICS = ['T-Count', 'T-Depth', 'Qubit-Count', 'Add Size', 'l']
        METRIC_NAME = METRICS[METRIC_INDEX]

        best_td_y = []
        best_tc_y = []
        best_qc_y = []
        best_qc_2_y = []
        # return TC, TD, QC, add_size, l

        plt.rc('font', family='Times New Roman', weight='bold', size=12)
        plt.rcParams['axes.spines.left'] = True
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.rcParams['axes.spines.bottom'] = True
        plt.figure(figsize=(4.5, 6.5))
        plt.xticks(np.arange(0, MAX_N + 1, 10))

        for N in range(3, max_N + 1):
            print(N)
            max_m = N
            Previous_Work = Our_Design(N, 3)
            print(f"Wallace_Mult_M_3 {Previous_Work}")
            Largest_M_Mul = Our_Design(N, max_m)
            print(f"Wallace_Mult_M_Largest {Largest_M_Mul}")
            Best_Rate_M_Mul = Our_Design(N, BestRate_M(N))
            print(f"Wallace_Mult_M_BestR {Best_Rate_M_Mul}")
            # Brute Force
            # 0,TC;      1,TD;       2,QC;   3,add_size; 4,max_l
            Best_M_Mul = Our_Design(N, Best_M(N, 1))
            print(f"Wallace_Mult_M_BestR {Best_M_Mul}")

            x.append(N)
            dp_td_y.append(dp(N))

            for i, add_size_fn_name in enumerate(add_size_fns):
                fn = add_size_fns[add_size_fn_name]
                Previous_Work_y[i].append(fn(Previous_Work)[METRIC_INDEX])
                Largest_M_Mul_y[i].append(fn(Largest_M_Mul)[METRIC_INDEX])
                Best_Rate_M_Mul_y[i].append(fn(Best_Rate_M_Mul)[METRIC_INDEX])
                Best_M_Mul_y[i].append(fn(Best_M_Mul)[METRIC_INDEX])
            print(Previous_Work_y)

            best_td_y.append(14 * N - 14)
            best_tc_y.append(18 * N * N - 24 * N)
            # best_qc_y.append(4 * N + 1)
            best_qc_2_y.append(2*N*N + 3*N)

        if METRICS[METRIC_INDEX] == 'T-Count':
            plt.plot(x, best_tc_y, label=f"Previous Work", color='grey', linestyle='--')

        if METRICS[METRIC_INDEX] == 'T-Depth':
            plt.plot(x, best_td_y, label=f"Previous Work", color='grey', linestyle='--')


        if METRICS[METRIC_INDEX] == 'Qubit-Count':
            # plt.plot(x, best_qc_y, label=f"Best Previous QC", color='red')
            plt.plot(x, best_qc_2_y, label=f"Previous Work", color='grey', linestyle='--')

        for y_name, y in [('Strategy 1', Largest_M_Mul_y),
                          ('Strategy 2', Best_Rate_M_Mul_y),
                          ('Strategy 3', Best_M_Mul_y)]:
        # print(Previous_Work_y)
        # for y_name, y in [('Previous Work', Previous_Work_y)]:
            for i, add_size_fn_name in enumerate(add_size_fns):

                if y_name == 'Strategy 1' and add_size_fn_name == 'Optimal TD':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-', color='mediumseagreen')
                elif y_name == 'Strategy 1' and add_size_fn_name == 'Optimal TC':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-', color='red')
                elif y_name == 'Strategy 2' and add_size_fn_name == 'Optimal TD':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-.', color='turquoise')
                elif y_name == 'Strategy 2' and add_size_fn_name == 'Optimal TC':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-.', color='dodgerblue')
                elif y_name == 'Strategy 3' and add_size_fn_name == 'Optimal TD':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='--', color='gold')
                elif y_name == 'Strategy 3' and add_size_fn_name == 'Optimal TC':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='--',   color='darkorange')
                else:
                    plt.plot(x, y[i], label= f"{y_name} ({add_size_fn_name})")

        # plt.legend()
        # leg = plt.legend(fontsize=11.5)
        #
        # leg.get_frame().set_alpha(0)
        plt.xlabel('N', fontsize=14, fontweight='bold')
        plt.ylabel(f'{METRICS[METRIC_INDEX]}', fontsize=13, fontweight='bold')
        # plt.title(f'Our Design - {METRICS[METRIC_INDEX]}', fontsize=16, fontweight='bold')
        plt.subplots_adjust(left=0.2, top=0.95, right=0.95)
        plt.savefig(f'Our Design {METRICS[METRIC_INDEX]}.pdf')
        plt.show()
        # if METRICS[METRIC_INDEX] == 'T-Count':
        #     input()

    # 110 addsize Cost TC= 101256 TD= 193 QC= 25514
    # This Cost TC= 101256 TD= 193 QC= 25514
