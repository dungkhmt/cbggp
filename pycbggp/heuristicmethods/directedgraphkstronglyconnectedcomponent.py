import random
from CBSGG import Graph,DirectedGraph

def generate_directed_graph_nb_strongly_connected_components(nb_nodes, nb_edges, nb_strongly_connected_components):
    scc = []
    m_max = 0

    # Find first position == x
    def Bsearch(l, r, x):
        while (l < r):
            mid = (l + r) // 2
            if (scc[mid] >= x):
                r = mid
            else:
                l = mid + 1
        return l

    # scc[i]--, scc[j]++
    def load(i, j, n):
        nonlocal m_max
        m_max -= scc[i] * (scc[i] - 1)
        m_max -= scc[j] * (scc[j] - 1)
        m_max -= scc[i] * (n - scc[i])
        m_max -= scc[j] * (n - scc[j])
        
        scc[i] -= 1
        scc[j] += 1
        
        m_max += scc[i] * (scc[i] - 1)
        m_max += scc[j] * (scc[j] - 1)
        m_max += scc[i] * (n - scc[i])
        m_max += scc[j] * (n - scc[j])

    # Gen k SCC
    def gen_k_scc(n, m, k):
        nonlocal m_max
        for i in range(k - 1, -1, -1):
            scc.append(n // k + (1 if i < (n % k) else 0))
            
        if (k == 1):
            return 
            
        tmp = 0
        for i in range(0, k):
            m_max += scc[i] * (scc[i] - 1)
            tmp += scc[i] * (n - scc[i])
            
        m_max += tmp // 2

        while (m_max < m - (k - 1)):
            l = Bsearch(0, k - 1, scc[k - 1])
            if (scc[l - 1] == 1):
                if (l == k - 1):
                    return 
                load(l, k - 1, n) 
            else:
                if (l - 2 < 0 or scc[l - 2] == 1):
                    load(l - 1, k - 1, n)  
                else :
                    l0 = Bsearch(0, l - 2, scc[l - 2])
                    load(l0, l - 1, n)  
        random.shuffle(scc)
    # Input
    n = nb_nodes
    m = nb_edges
    k = nb_strongly_connected_components

    gen_k_scc(n, m, k)

    edge = []
    if (m < n and not(m == n - 1 and k == n)) :
        print(-1)
    else :
        # Tạo khung cho k scc
        # for i in range(0, k):
        #     print(scc[i])
        # print()
        m0 = 0
        l = [0] * k
        r = [0] * k
        cnt = 0
        mp = set()
        edge_scc = [[] for _ in range(0, k)]
        # Tạo sẵn 1 chu trình cho mỗi scc
        for i in range(0, k):
            l[i] = cnt + 1;
            r[i] = cnt + scc[i]
            cnt += scc[i]
            if (scc[i] == 1):
                continue
            for j in range(l[i], r[i] + 1):
                u = j
                v = j + 1 if (j < r[i]) else l[i]
                mp.add((u, v))
                edge.append((u, v))
                m0 += 1
            # Với những cụm scc <= 100 đỉnh thì tạo danh sách cạnh ứng viên 
            # Với những cum scc > 100 thì sẽ tạo cạnh random
            if (scc[i] <= 100):
                for u in range(l[i], r[i] + 1):
                    for v in range(l[i], r[i] + 1):
                        if (u != v and u + 1 != v):
                            if (u == r[i] and v == l[i]):
                                continue
                            edge_scc[i].append((u, v))
                random.shuffle(edge_scc[i])
        # Tạo trước 1 khung để kết nối k_scc, đảm bảo đồ thị liên thông yếu
        for i in range(1, k):
            j = random.randint(0, i - 1)
            u = random.randint(l[j], r[j])
            v = random.randint(l[i], r[i])
            mp.add((u, v))
            edge.append((u, v))
            m0 += 1
        # Tạo các cạnh nối các thành phần scc với nhau
        # Chỉ nối chỉ số bé tới lớn để đảm bảo không tạo ra chu trình mới ngoài k scc đã có
        # Nếu n <= 3000 thì tạo danh sách cách ứng viên nối giữa các scc 
        # Nếu n > 3000 thì sẽ random chọn cạnh

        edge_dag = []
        if (n <= 3000):
            for i in range(0, k):
                for u in range(l[i], r[i] + 1):
                    for v in range(r[i] + 1, n + 1):
                        if ((u, v) not in mp):
                            edge_dag.append((u, v))
        
        # Thêm cạnh đồ thị vào để đủ m cạnh
        rd_scc = list(range(0, k))
        random.shuffle(rd_scc)
        while (m0 < m):
            t = random.randint(0, 1)
            # t = 0 thì thêm cạnh nội bộ của 1 scc bất kì
            # t = 1 thì thêm cạnh giữa các scc
            if (t == 0):
                if (len(rd_scc) == 0):
                    continue
                id = random.randint(0, len(rd_scc) - 1)
                if (scc[rd_scc[id]] <= 100):
                    if (len(edge_scc[rd_scc[id]]) > 0):
                        u, v = edge_scc[rd_scc[id]].pop()
                        edge.append((u, v))
                        m0 += 1
                    else :
                        last_val = rd_scc.pop()
                        if id < len(rd_scc):
                            rd_scc[id] = last_val
                else:
                    idx = rd_scc[id]
                    u = random.randint(l[idx], r[idx])
                    v = random.randint(l[idx], r[idx])
                    if (u != v and (u, v) not in mp):
                        m0 += 1;
                        edge.append((u, v))
                        mp.add((u, v))
            else :
                if (n <= 3000):
                    if (len(edge_dag) > 0):
                        u, v = edge_dag.pop()
                        edge.append((u, v))
                        m0 += 1
                elif (k > 1):
                    i = random.randint(0, k - 2)
                    j = random.randint(i + 1, k - 1)
                    u = random.randint(l[i], r[i])
                    v = random.randint(l[j], r[j])
                    if ((u, v) not in mp):
                        m0 += 1;
                        edge.append((u, v))
                        mp.add((u, v))
    sf_idx = list(range(1, n + 1))
    random.shuffle(sf_idx)
    random.shuffle(edge)

    Res = DirectedGraph(n)
    for (u, v) in edge:
        Res.AddEdge(sf_idx[u - 1] - 1, sf_idx[v - 1] - 1)
    return Res
