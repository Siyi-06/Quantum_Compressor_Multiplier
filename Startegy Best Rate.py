import math

# Parameters
N = 8  # Size of the multiplier (N x N)
max_m = N  # Maximum value of m for m:k compressors
max_layers = N  # Maximum number of layers in the Wallace tree, not include the last Carry Propagation Adder
max_columns = 2*N-1 #Maximum columns in Wallace Tree
P_value = [[0 for l in range(max_layers)]for j in range(max_columns)]
#print(P_value)
def TIME(m):
    time=math.floor(math.log(m_largest, 3)) + 2 * math.floor(math.log(m_largest, 2)) - 1
    return time

# First Layer: Initialize the partial product count
for j in range(N):
    P_value[j][0] = j + 1
    # print("!", j, j + 1)
    if j!=N-1:
        P_value[max_columns - 1 - j][0] = j + 1
        # print("!",max_columns -1- j, j + 1)
print(P_value)

#Cost of First Layer
TC=4*N*N
TD=N
QC=2*N+N*N
print("Cost of First Layer:")
print("TC=",TC,"TD=",TD,"QC=",QC)

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
    print("flag!!!",flag_value)
    if flag_value==0:
        print("max_l=",l)#maximum level 3 layers, not include CPA
        print("add_size=", add_size)
        break
    for j in range(max_columns):
        # Apply compressors and update the partial products
        if P_value[j][l]>1:
            #m=P_value[j][l]#use the largest compressor always
            k = math.floor(math.log2(P_value[j][l])) + 1
            m =math.floor(math.pow( 2, k )-1)
            #print("P_value[",j,",",l,"]",P_value[j][l],"m",m,"k",k)
            if m>P_value[j][l]:
                k=k-1
                m=math.floor(math.pow( 2, k )-1)
            if m==2:
                k=2
            if P_value[j][l]==2:
                m=2
                k=2
            if m<2:
                print("no compressor here!","m",m)
                continue
            print("P_value[", j, ",", l, "]", P_value[j][l], "m", m, "k", k)

            ###Update Cost
            if m==0:
                continue
            if m_largest<m:
                m_largest=m
            if m>k:
                TC = TC + 4 * (m - k)
                QC = QC + (m - k)
            if m==2:
                TC=TC+4
                QC=QC+1
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
                        print("carry propagation")
                        print("P_value[", j+ii, ",", l+1 , "]", P_value[j+ii][l+1])
    #Update Cost
    print("Previous TD", TD)
    TD=TD+math.floor(math.log(m_largest,3))+2*math.floor(math.log(m_largest,2))-2
    print("m_largest", m_largest)
    print("TD",TD)
print(P_value)
print("Cost", "TC=",TC,"TD=",TD,"QC=",QC)

