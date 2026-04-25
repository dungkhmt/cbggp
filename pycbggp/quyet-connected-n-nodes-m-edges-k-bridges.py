import sys
import random
from CBSGG import Graph, Edge
def delta_max(x, y, t):
    old_x = x * (x - 1) // 2
    old_y = y * (y - 1) // 2
    x -= t
    y += t
    new_x = x * (x - 1) // 2
    new_y = y * (y - 1) // 2
    return new_x + new_y - old_x - old_y
def delta_min(x, y, t):
    old_x = x if x >= 3 else 0
    old_y = y if y >= 3 else 0
    x -= t;
    y += t;
    new_x = x if x >= 3 else 0
    new_y = y if y >= 3 else 0
    return new_x + new_y - old_x - old_y
    
def partition_vertices(n, m, k):
    c = k + 1
    if k > n - 1 or k == n - 2:
        return None
    if m < n - 1 or m > k + (n - k) * (n - k - 1) // 2:
        return None
    s = [1] * c
    s[0] += n - c
    max_edges = 0
    min_edges = 0
    for sz in s:
        max_edges += sz * (sz - 1) // 2
        if (sz >= 3):
            min_edges += sz
    in_PQ = [0] * (c)
    in_PQ[0] = 1
    PQ = [0] if s[0] > 2 else []
    for step in range(100000): 
        u = random.choice(PQ)
        v = random.randint(0, c - 1)
        if (u == v):
            continue
        cnt_move = 1
        if (s[v] == 1):
            cnt_move = 2
        if (s[u] - cnt_move >= 3 or s[u] - cnt_move == 1):
            delta_max_edges = delta_max(s[u], s[v], cnt_move)
            delta_min_edges = delta_min(s[u], s[v], cnt_move)
            if (max_edges + delta_max_edges >= m - k and min_edges + delta_min_edges <= m - k):
                max_edges += delta_max_edges
                min_edges += delta_min_edges
                s[u] -= cnt_move
                s[v] += cnt_move
                if (s[v] >= 3 and in_PQ[v] == 0):
                    in_PQ[v] = 1;
                    PQ.append(v)
    return s
n, m, k = map(int, input().split())

S = partition_vertices(n, m, k)

if S is None : 
    print(-1)
else:
    edge = []
    L = [0] * (k + 1)
    R = [0] * (k + 1)
    cnt = 0
    for i in range(0, k + 1):
        L[i] = cnt;
        R[i] = cnt + S[i] - 1
        cnt += S[i]
    # Bridges
    for i in range(1, k + 1):
        u = random.randint(L[i], R[i])
        j = random.randint(0, i - 1)
        v = random.randint(L[j], R[j])
        edge.append((u, v))
    ec = [0] * (k + 1)
    ST = []
    num_e = m - k
    for i in range(0, k + 1):
        if (S[i] >= 3):
            ec[i] = S[i]
            ST.append(i)
            num_e -= S[i]
    while (num_e > 0 and len(ST) > 0):
        sz = len(ST)
        i = random.randint(0, sz - 1)
        j = ST[i]
        if (ec[j] == S[j] * (S[j] - 1) // 2):
            tmp = ST[i]
            ST[i] = ST[sz - 1]
            ST[sz - 1] = tmp
            ST.pop()
            continue
        ec[j] += 1
        num_e -= 1
    # generate_ecc_edges
    for i in range(0, k + 1):
        if (S[i] == 1):
            continue
        V = list(range(L[i], R[i] + 1))
        random.shuffle(V)
        sz = S[i]
        mp = set()
        for j in range(0, sz):
            u = V[j]
            v = V[(j + 1) % sz]
            edge.append((u, v))
            if (u > v):
                u, v = v, u
            mp.add((u, v))
        need_edge = ec[i] - sz
        max_edges = sz * (sz - 1) // 2
        if (need_edge < max_edges // 2):
            while (need_edge > 0):
                u = random.choice(V)
                v = random.choice(V)
                if (u == v):
                    continue
                if (u > v):
                    u, v = v, u
                if ((u, v) not in mp):
                    mp.add((u, v))
                    sw = random.randint(0, 1)
                    if (sw == 0):
                        u, v = v, u
                    edge.append((u, v))
                    need_edge -= 1
        else :
            all_edge = []
            for i in range(0, sz):
                for j in range(i + 1, sz):
                    u = V[i]
                    v = V[j]
                    if (u > v):
                        u, v = v, u
                    if ((u, v) not in mp):
                        sw = random.randint(0, 1)
                        if (sw == 0):
                            u, v = v, u
                        all_edge.append((u, v))
            random.shuffle(all_edge)
            for i in range(0, need_edge):
                u, v = all_edge[i]
                edge.append((u, v))
    # Random vertex
    Idx = list(range(1, n + 1))
    random.shuffle(Idx)
    random.shuffle(edge)
    # print(n, end = " ")
    # print(m)
    # for (u, v) in edge:
    #     print(Idx[u], end = " ")
    #     print(Idx[v])
    Res = Graph(n)
    for (u, v) in edge:
        Res.AddEdge(Idx[u] - 1, Idx[v] - 1)

    Res.Print()



