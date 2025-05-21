import math


def get_k(m):
    k = math.floor(math.log2(m)) + 1
    return k


def get_cost(m):
    if m == 2 or m == 3:
        return 1
    else:
        return math.ceil(math.log(m, 3)) + 2 * math.floor(math.log2(m)) - 2


def recursive(N, dp, dp_path):
    if N <= 2 or dp[N] != math.inf:
        return dp[N]
    else:
        for m in range(N, 2, -1):
            k = get_k(m)
            diff = m - k
            next_N = N - diff
            cur_cost = get_cost(m)

            return_cost = recursive(next_N, dp, dp_path)
            if return_cost != math.inf:
                new_cost = return_cost + cur_cost
                if new_cost < dp[N]:
                    dp_path[N] = m
                    dp[N] = new_cost
        return dp[N]


def dp(N, show_print=False, cpa = 'optimal_td'):
    dp = [math.inf] * (N + 1)
    dp_m = [-1] * (N + 1)
    dp[2] = 0

    recursive(N, dp, dp_m)

    if cpa == 'optimal_td':
        cpa = math.ceil(math.log2(N) + 1)
    elif cpa == 'optimal_tc':
        cpa = N

    total_td = dp[N] + cpa + N
    
    if show_print:
        print("TD:", dp[N], cpa, dp[N] + cpa + N)
        print_optimal_path(N, dp_m)

    return total_td

def print_optimal_path(N, dp_m):
    cur = N

    while cur > 2:
        m = dp_m[cur]
        k = get_k(m)
        diff = m - k
        print(cur, m, get_cost(m))
        cur -= diff
    print(cur, m, get_cost(m))


if __name__ == "__main__":
    dp(64, show_print=True)
