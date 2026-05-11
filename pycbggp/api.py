from CBSGG import Graph,DirectedGraph
import random
# undirected graphs
def gen_undirected_graph(nb_nodes, nb_edges, nb_connected_components, nb_bridges, nb_articulation_points):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  
 
def gen_undirected_connected_graph_no_bridge(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph_no_articulation_point(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_connected_graph_no_bridge_no_articulation_point(nb_nodes, nb_edges):
 G = Graph(nb_nodes)
 return G  
def gen_undirected_connected_graph_nb_bridges(nb_nodes, nb_edges, nb_bridges):
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
    n = nb_nodes
    m = nb_edges
    k = nb_bridges

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
        return Res
# special undirected graphs 
def gen_undirected_complete_graph(nb_nodes):
 G = Graph(nb_nodes)
 return G  

def gen_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 G = Graph(nb_nodes)
 return G  

def gen_connected_bipartite_graph(nb_left_nodes, nb_right_nodes,nb_edges):
 G = Graph(nb_nodes)
 return G  


# undirected trees
def gen_undirected_tree(nb_nodes):
 G = Graph(nb_nodes)
 return G  

def gen_undirected_tree_bounded_diameter_degree(nb_nodes, ub_deg, ub_diameter):
 G = Graph(nb_nodes)
 return G  
  
# directed graphs 
def gen_directed_graph(nb_nodes, nb_edges):
 G = DirectedGraph(nb_nodes)
 return G  

def gen_directed_strongly_connected_graph(nb_nodes, nb_edges):
 G = DirectedGraph(nb_nodes)
 return G 
 
def gen_directed_graph_nb_strongly_connected_components(nb_nodes, nb_edges, nb_strongly_connected_components):
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
 

# main    
G = gen_undirected_connected_graph_nb_bridges(7, 8, 2)
G.Print()

 
 