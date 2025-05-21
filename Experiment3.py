import matplotlib.pyplot as plt
import math
import numpy as np
from Startegy_Dynamic_Programming import dp

def Highest_Strategy(N):
    # Parameters
    #N = 1024  # Size of the multiplier (N x N)
    max_m = N  # Maximum value of m for m:k compressors
    max_layers = N  # Maximum number of layers in the Wallace tree, not include the last Carry Propagation Adder
    max_columns = 2*N-1 #Maximum columns in Wallace Tree
    P_value = [[0 for l in range(max_layers)]for j in range(max_columns)]
    #print(P_value)

    # First Layer: Initialize the partial product count
    for j in range(N):
        P_value[j][0] = j + 1
        # print("!", j, j + 1)
        if j!=N-1:
            P_value[max_columns - 1 - j][0] = j + 1
            # print("!",max_columns -1- j, j + 1)
    #print(P_value)

    #Cost of First Layer
    TC=4*N*N
    TD=N
    QC=2*N+N*N
    #print("Cost of First Layer:")
    #print("TC=",TC,"TD=",TD,"QC=",QC)

    # Compressor constraints and partial product updates
    flag_value=3
    add_size=100
    m_largest=100
    for l in range(max_layers):
        #flag for check last layer
        m_largest=0
        flag_value=0
        add_size=0
        for j in range(max_columns):
            if add_size==0:
                if P_value[j][l]==2:
                    add_size=max_columns-j
            if P_value[j][l] > 2:
                flag_value = 1
        #print("flag!!!",flag_value)
        if flag_value==0:
            #print("max_l=",l)#maximum level 3 layers, not include CPA
            #print("add_size=", add_size)
            break
        for j in range(max_columns):
            # Apply compressors and update the partial products
            if P_value[j][l]>1:
                m=P_value[j][l]#use the largest compressor always
                k = math.floor(math.log2(m)) + 1
                #print("P_value[",j,",",l,"]",P_value[j][l],"m",m,"k",k)
                ###Update Cost
                if m_largest<m:
                    m_largest=m
                TC = TC + 4 * (m - math.floor(k/2))
                QC = QC + (m- math.floor(k/2))

                # Update the partial product count for column j after applying the compressor
                # Reduce the partial products at column j
                if l < (max_layers-1):
                    # Reduce the partial products at column j
                    P_value[j][l + 1] = P_value[j][l] - m + 1 + P_value[j][l + 1]
                    # print("Compression")
                    # print("P_value[", j, ",", l, "]", P_value[j][l])
                    # print("P_value[", j, ",", l+1, "]", P_value[j][l + 1])

                    # Propagate carry-out bits to the next columns (j+1 to j+k-1)
                    for ii in range(1,k):
                        if j + ii < 2 * N - 1:
                            # print("Before carry propagation")
                            # print("P_value[", j + ii, ",", l + 1, "]", P_value[j + ii][l + 1])
                            P_value[j+ii][l+1] = P_value[j+ii][l+1]+1
                            #print("carry propagation")
                            # print("P_value[", j+ii, ",", l+1 , "]", P_value[j+ii][l+1])
        #Update Cost
        #print("Previous TD", TD)
        TD=TD+math.ceil(math.log(m_largest,3))+2*math.floor(math.log(m_largest,2))-2
        #print("m_largest", m_largest)
        #print("TD",TD)
    #print(P_value)
    #print("Cost", "TC=",TC,"TD=",TD,"QC=",QC)
    return TC, TD, QC, add_size, l


def Best_Rate_Strategy(N):
    # Parameters
    #N = 1024  # Size of the multiplier (N x N)
    max_m = N  # Maximum value of m for m:k compressors
    max_layers = N  # Maximum number of layers in the Wallace tree, not include the last Carry Propagation Adder
    max_columns = 2 * N - 1  # Maximum columns in Wallace Tree
    P_value = [[0 for l in range(max_layers)] for j in range(max_columns)]

    # print(P_value)
    # def TIME(m):
    #     time = math.ceil(math.log(m_largest, 3)) + 2 * math.floor(math.log(m_largest, 2)) - 1
    #     return time

    # First Layer: Initialize the partial product count
    for j in range(N):
        P_value[j][0] = j + 1
        # print("!", j, j + 1)
        if j != N - 1:
            P_value[max_columns - 1 - j][0] = j + 1
            # print("!",max_columns -1- j, j + 1)
    #print(P_value)

    # Cost of First Layer
    TC = 4 * N * N
    TD = N
    QC = 2 * N + N * N
    #print("Cost of First Layer:")
    #print("TC=", TC, "TD=", TD, "QC=", QC)

    # Compressor constraints and partial product updates
    flag_value = 3
    add_size = 100
    m_largest = 100
    for l in range(max_layers):
        # flag for check last layer
        m_largest = 0
        flag_value = 0
        add_size = 0
        for j in range(max_columns):
            if add_size == 0:
                if P_value[j][l] == 2:
                    add_size = max_columns - j
            if P_value[j][l] > 2:
                flag_value = 1
        #print("flag!!!", flag_value)
        if flag_value == 0:
            #print("max_l=", l)  # maximum level 3 layers, not include CPA
            #print("add_size=", add_size)
            break
        for j in range(max_columns):
            # Apply compressors and update the partial products
            if P_value[j][l] > 1:
                # m=P_value[j][l]#use the largest compressor always
                k = math.floor(math.log2(P_value[j][l])) + 1
                m = math.floor(math.pow(2, k) - 1)
                # print("P_value[",j,",",l,"]",P_value[j][l],"m",m,"k",k)
                if m > P_value[j][l]:
                    k = k - 1
                    m = math.floor(math.pow(2, k) - 1)
                if m == 2:
                    k = 2
                if P_value[j][l] == 2:
                    m = 2
                    k = 2
                if m < 2:
                    #print("no compressor here!", "m", m)
                    continue
                #print("P_value[", j, ",", l, "]", P_value[j][l], "m", m, "k", k)

                ###Update Cost
                if m == 0:
                    continue
                if m_largest < m:
                    m_largest = m
                if m > k:
                    TC = TC + 4 * (m - k)
                    QC = QC + (m - k)
                if m == 2:
                    TC = TC + 4
                    QC = QC + 1
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
                            #print("carry propagation")
                            #print("P_value[", j + ii, ",", l + 1, "]", P_value[j + ii][l + 1])
        # Update Cost
        #print("Previous TD", TD)
        TD = TD + math.ceil(math.log(m_largest,3)) + 2 * math.floor(math.log(m_largest, 2)) - 2
        #print("m_largest", m_largest)
        #print("TD", TD)
    #print(P_value)
    #print("Cost", "TC=", TC, "TD=", TD, "QC=", QC)
    return TC, TD, QC, add_size, l

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

def DP(N):
    max_m = N  # 最大压缩器的输入数
    partial_products = [i for i in range(1, N + 1)]
    for i in range(1, N):
        partial_products.append(N - i)
    #print(partial_products)
    max_columns = len(partial_products)
    max_layers = math.floor((len(partial_products) + 1) / 2)
    dp_TD = [math.inf] * max_layers
    Compressor_depth = [0] * (max_layers + 1)
    Compressor_depth[0] = 0
    Compressor_depth[1] = 0
    Compressor_depth[2] = 1
    Compressor_depth[3] = 1
    for i in range(4, max_layers + 1):
        #print(i)
        Compressor_depth[i] = math.ceil(math.log(i, 3)) + 2 * math.floor(math.log(i, 2)) - 2
    #print(Compressor_depth)
    # Cost of First Layer
    TC = 4 * N * N
    TD = N
    QC = 2 * N + N * N

    # dp_TD[0]=TD
    for i in range(1, max_layers + 1):
        flag = 0
        #print("partial_products", partial_products)
        max_p = max(partial_products)
        #print("max_p", max_p)
        if max_p - 3 * i + i <= 2:
            # print("!")
            flag = 1
            continue
        #print(i)
        #print(flag)
        #print(max_p - 3 * i + i)
        if flag == 0:
            dp_TD[i] = N + min(i * Compressor_depth[3] + Compressor_depth[max_p - 3 * i + i] + Compressor_depth[
                math.floor(math.log2(max_p - 3 * i + i))],
                               Compressor_depth[max_p] + Compressor_depth[math.floor(math.log2(max_p))])
        else:
            #print("flag")
            dp_TD[i] = N + Compressor_depth[max_p] + Compressor_depth[math.floor(math.log2(max_p))]
            #print()
        #print(N, i * Compressor_depth[3], Compressor_depth[max_p - 3 * i + i])
        #print(N + i * Compressor_depth[3] + Compressor_depth[max_p - 3 * i + i])
        #print(N + Compressor_depth[max_p] + Compressor_depth[math.floor(math.log2(max_p))])
        #print("dp_TD", dp_TD)
        if max_p - 3 * i + i < 3:
            break
        # 要么分解，要么一步
    #print("dp_TD", dp_TD)
    return (min(dp_TD))

def Our_Design(N,max_m):
    # Parameters
    #N = 100  # Size of the multiplier (N x N)
    #max_m = 15  # Maximum value of m for m:k compressors
    max_layers = N  # Maximum number of layers in the Wallace tree, not include the last Carry Propagation Adder
    max_columns = 2*N-1 #Maximum columns in Wallace Tree
    P_value = [[0 for l in range(max_layers)]for j in range(max_columns)]
    #print(P_value)

    # First Layer: Initialize the partial product count
    for j in range(N):
        P_value[j][0] = j + 1
        # print("!", j, j + 1)
        if j!=N-1:
            P_value[max_columns - 1 - j][0] = j + 1
            # print("!",max_columns -1- j, j + 1)
    #print(P_value)

    #Cost of First Layer
    TC=4*N*N
    TD=N
    QC=2*N+N*N
    #print("Cost of First Layer:")
    #print("TC=",TC,"TD=",TD,"QC=",QC)

    # Compressor constraints and partial product updates
    flag_value=3
    add_size=100
    m_largest=100
    m=10000
    k=10000
    for l in range(max_layers):
        #flag for check last layer
        m_largest=0
        flag_value=0
        add_size=0
        for j in range(max_columns):
            if add_size==0:
                if P_value[j][l]==2:
                    add_size=max_columns-j
                    break
        for j in range(max_columns):
            if P_value[j][l] > 2:
                flag_value = 1
        #print("flag!!!",flag_value)
        if flag_value==0:
            #print("max_l=",l)#maximum level 3 layers, not include CPA
            #print("add_size=", add_size)
            break
        for j in range(max_columns):
            # Apply compressors and update the partial products
            if P_value[j][l] > 1:
                m = P_value[j][l]  # use the largest compressor always
                k = math.floor(math.log2(m)) + 1
            if P_value[j][l]<2:
                m=0
                k=0
                #print("P_value[", j, ",", l, "]", P_value[j][l], "m", m, "k", k)
                #print("No compressor")
                continue
            if l == 0:
                if P_value[j][l] > max_m:
                    m = max_m
                    k = math.floor(math.log2(m)) + 1
                    #print("YES")
                #print("P_value[", j, ",", l, "]", P_value[j][l], "m", m, "k", k)
            if l>0:
                if P_value[j][l] > 2:
                    m=3
                    k=2
                if P_value[j][l]==2:
                    m=2
                    k=2
                if P_value[j][l] < 2:
                    #print("No Compressor Here!")
                    m=0
                    k=0
                    continue
                #print("P_value[",j,",",l,"]",P_value[j][l],"m",m,"k",k)
                ###Update Cost
            #print("m_largest", m_largest)
            if m_largest<m:
                m_largest=m
                #print("m_largest",m_largest)
            if m > k:
                TC = TC + 4 * (m - k)
                QC = QC + (m - k)
            if m == k:
                TC = TC + 4
                QC = QC + 1
                #print("!!!B2")

                # Update the partial product count for column j after applying the compressor
                # Reduce the partial products at column j
            if l < (max_layers-1):
                    # Reduce the partial products at column j
                P_value[j][l + 1] = P_value[j][l] - m + 1 + P_value[j][l + 1]
                    # print("Compression")
                    # print("P_value[", j, ",", l, "]", P_value[j][l])
                    # print("P_value[", j, ",", l+1, "]", P_value[j][l + 1])

                    # Propagate carry-out bits to the next columns (j+1 to j+k-1)
                for ii in range(1,k):
                    if j + ii < 2 * N - 1:
                            # print("Before carry propagation")
                            # print("P_value[", j + ii, ",", l + 1, "]", P_value[j + ii][l + 1])
                        P_value[j+ii][l+1] = P_value[j+ii][l+1]+1
                        #print("carry propagation")
                            # print("P_value[", j+ii, ",", l+1 , "]", P_value[j+ii][l+1])
        #Update Cost
        #print("Previous TD", TD)
        TD=TD+math.ceil(math.log(m_largest,3))+2*math.floor(math.log(m_largest,2))-2
        #print("m_largest", m_largest)
        #print("TD",TD)
    #print(P_value)
    #print("Cost", "TC=",TC,"TD=",TD,"QC=",QC)
    #TC,TD,QC,add_size,max_l
    return TC,TD,QC,add_size,l

def Best_M(N,metric):
    flag=[math.inf,math.inf,math.inf,math.inf,math.inf]
    Best_M=0
    for M in range(3,N+1):
        #print(N,M)
        value=Our_Design(N, M)
        if value[metric]<flag[metric]:
            flag[metric]=value[metric]
            Best_M=M
    #print('Best_M',Best_M)
    return Best_M


if __name__ == "__main__":
    add_size_fns = {'Optimal TD': optimal_td,
                    'Optimal TC': gidney_tc}
    max_N=39

    for METRIC_INDEX in [0]:

        best_td_y = []
        best_tc_y = []
        best_qc_y = []

        td_dp_optimal_td_y = []
        td_dp_optimal_tc_y = []

        x = []
        Largest_M_Mul_y = [[] for _ in range(len(add_size_fns))]
        Best_Rate_M_Mul_y = [[] for _ in range(len(add_size_fns))]
        Best_M_Mul_y = [[] for _ in range(len(add_size_fns))]
        Best_Depth_y = []

        plt.rc('font', family='Times New Roman', weight='bold', size=12)
        plt.rcParams['axes.spines.left'] = True
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
        plt.rcParams['axes.spines.bottom'] = True
        plt.figure(figsize=(4.5, 6.5))
        plt.xticks(np.arange(0, max_N + 1, 10))

        Previous_Work_y = [[] for _ in range(len(add_size_fns))]


        METRICS = ['T-Count', 'T-Depth', 'Qubit-Count', 'Add Size', 'l']

        # 计算最优 T-depth
        for N in range(3,max_N+1):
            print(N)
            max_m = N
            x.append(N)
            td_dp_optimal_td_y.append(dp(N))
            td_dp_optimal_tc_y.append(dp(N, cpa='optimal_tc'))
            #Previous_Work = Our_Design(N, 3)
            #print(f"Wallace_Mult_M_3 {Previous_Work}")
            Largest_M_Mul = Highest_Strategy(N)
            Largest_M_Mul_y.append(Largest_M_Mul[METRIC_INDEX])
            print(f"Wallace_Mult_M_Largest {Largest_M_Mul}")

            Best_Rate_M_Mul=Best_Rate_Strategy(N)
            Best_Rate_M_Mul_y.append(Best_Rate_M_Mul[METRIC_INDEX])
            print(f"Wallace_Mult_M_BestR {Best_Rate_M_Mul}")
            # Brute Force
            # 0,TC;      1,TD;       2,QC;   3,add_size; 4,max_l
            Best_M_Mul = Our_Design(N, Best_M(N, 1))
            Best_M_Mul_y.append(Best_M_Mul[METRIC_INDEX])
            print(f"Wallace_Mult_OurDesign {Best_M_Mul}")
            #Brute Force
            #0,TC;      1,TD;       2,QC;   3,add_size; 4,max_l
            if N>4:
                best_depth = DP(N)
                print(f"Best {best_depth}")
                Best_Depth_y.append(best_depth)
            else:
                Best_Depth_y.append(0)

            for i, add_size_fn_name in enumerate(add_size_fns):
                fn = add_size_fns[add_size_fn_name]
                Largest_M_Mul_y[i].append(fn(Largest_M_Mul)[METRIC_INDEX])
                Best_Rate_M_Mul_y[i].append(fn(Best_Rate_M_Mul)[METRIC_INDEX])
                Best_M_Mul_y[i].append(fn(Best_M_Mul)[METRIC_INDEX])



            # best_td_y.append(13 * N - 13)
            # best_tc_y.append(18 * N * N - 24 * N)
            # best_qc_y.append(2 * N * N + N + 2 * N)

        # if METRICS[METRIC_INDEX] == 'T-Count':
        #     # plt.plot(x, best_tc_y, label=f"Best Previous TC", color='red')
        #
        # if METRICS[METRIC_INDEX] == 'T-Depth':
        #     # plt.plot(x, best_td_y, label=f"Best Previous TD", color='red')
        #
        # if METRICS[METRIC_INDEX] == 'Qubit-Count':
            # plt.plot(x, best_qc_y, label=f"Best Previous QC", color='red')
            # plt.plot(x, best_qc_y, label=f"Best Previous QC 2", color='red')

        for y_name, y in [('Construction 1', Best_M_Mul_y),
                          ('Strategy A', Largest_M_Mul_y),
                          ('Strategy B', Best_Rate_M_Mul_y)]:
        # print(Previous_Work_y)
        # for y_name, y in [('Previous Work', Previous_Work_y)]:
            for i, add_size_fn_name in enumerate(add_size_fns):
                # plt.plot(x, y[i], label= f"{y_name} ({add_size_fn_name})")
                if y_name == 'Strategy A' and add_size_fn_name == 'Optimal TD':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-', color='mediumseagreen')
                elif y_name == 'Strategy A' and add_size_fn_name == 'Optimal TC':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-', color='red')
                elif y_name == 'Strategy B' and add_size_fn_name == 'Optimal TD':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-.', color='turquoise')
                elif y_name == 'Strategy B' and add_size_fn_name == 'Optimal TC':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='-.', color='dodgerblue')
                elif y_name == 'Construction 1' and add_size_fn_name == 'Optimal TD':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='dotted', color='grey')
                elif y_name == 'Construction 1' and add_size_fn_name == 'Optimal TC':
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})", linestyle='--', color='grey')
                else:
                    plt.plot(x, y[i], label=f"{y_name} ({add_size_fn_name})")

        # plt.plot(x, td_dp_optimal_td_y, label=f"Strategy C (Optimal TD)", linestyle='-', color='gold')
        # plt.plot(x, td_dp_optimal_tc_y, label=f"Strategy C (Optimal TC)", linestyle='-', color='darkorange')

        plt.xlim(0,max_N)


        # if METRIC_INDEX == 1:
        #     plt.plot(x, Best_Depth_y, label=f"Best {METRICS[METRIC_INDEX]}")

        # plt.legend()
        # leg = plt.legend(fontsize=11.5)
        # leg.get_frame().set_alpha(0)

        plt.xlabel('N', fontsize=14, fontweight='bold')
        plt.ylabel(f'{METRICS[METRIC_INDEX]}', fontsize=13, fontweight='bold')

        plt.subplots_adjust(left=0.2, top=0.95, right=0.95)
        plt.savefig(f'More Compressors {METRICS[METRIC_INDEX]}.pdf')
        # plt.title(f'More Compressors - {METRICS[METRIC_INDEX]}')

        plt.show()