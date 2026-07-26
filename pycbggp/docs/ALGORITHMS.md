# Tài Liệu Thuật Toán Sinh Đồ Thị 

**Tác giả:** Nguyen Ngoc Tuan Anh  

---

## Đồ thị vô hướng liên thông
**File:** `pycbggp/constructivemethods/undirected_connected_graph.py`  
**Hàm:** `generate_undirected_connected_graph(nb_nodes, nb_edges)`

### Điều kiện khả thi
- `V - 1 ≤ E ≤ V*(V-1)/2`

### Thuật toán — Cây khung ngẫu nhiên + Bổ sung cạnh
1. **Cây khung:** Trộn ngẫu nhiên tất cả `V*(V-1)/2` cạnh ứng viên, duyệt với DSU (Disjoint Set Union). Chấp nhận `V-1` cạnh đầu tiên hợp nhất hai thành phần khác nhau → đảm bảo cây khung liên thông.
2. **Cạnh thêm:** Các cạnh còn lại sau bước DSU đã ở thứ tự ngẫu nhiên. Lấy `E-(V-1)` cạnh đầu tiên.
3. **Hoán vị nhãn:** Áp dụng hoán vị ngẫu nhiên cho toàn bộ ID đỉnh.

**Độ phức tạp:** O(V² + E·α(V)), trong đó α là hàm nghịch đảo Ackermann từ DSU.

---

## Đồ thị vô hướng liên thông không có cầu (2-Liên thông cạnh)
**File:** `pycbggp/constructivemethods/undirected_connected_graph_no_bridge.py`  
**Hàm:** `generate_undirected_connected_graph_no_bridge(nb_nodes, nb_edges)`

### Điều kiện khả thi
- `V ≥ 3`, `V ≤ E ≤ V*(V-1)/2`

### Thuật toán — Chu trình Hamilton + Cạnh thêm
1. **Chu trình Hamilton:** Trộn ngẫu nhiên tất cả đỉnh thành hoán vị `perm`. Thêm các cạnh `perm[0]-perm[1]`, `perm[1]-perm[2]`, ..., `perm[V-1]-perm[0]`. Tạo ra chu trình độ dài V gồm V cạnh — mỗi cạnh đều nằm trên chu trình nên không có cầu.
2. **Cạnh thêm:** Thu thập tất cả cặp đỉnh chưa có cạnh, trộn ngẫu nhiên, lấy `E-V` cạnh đầu tiên.

**Tính đúng đắn:** Thêm bất kỳ cạnh nào vào đồ thị 2-liên thông cạnh vẫn giữ nguyên tính 2-liên thông cạnh (chỉ tạo thêm chu trình, không thể tạo cầu mới).  
**Độ phức tạp:** O(V²)

---

## Đồ thị 2-liên thông đỉnh (Không có điểm khớp, Không có cầu)
**File:** `pycbggp/constructivemethods/biconnected_graph.py`  
**Hàm:** `generate_biconnected_graph(nb_nodes, nb_edges)`

Dùng cho cả hai bài toán:
- `gen_undirected_connected_graph_no_articulation_point` (Task 3)
- `gen_undirected_connected_graph_no_bridge_no_articulation_point` (Task 4)

Hai bài toán tương đương vì: **2-liên thông đỉnh ⟺ không có điểm khớp ⟹ không có cầu**.

### Điều kiện khả thi
- `V ≥ 3`, `V ≤ E ≤ V*(V-1)/2`

### Thuật toán — Phân tích tai mở (Open Ear Decomposition — Whitney, 1932)
Đồ thị là 2-liên thông đỉnh khi và chỉ khi nó có phân tích tai mở. Một **tai mở** là đường đi `u → w₁ → w₂ → ... → v` trong đó `u`, `v` đã có trong đồ thị (`u ≠ v`) và tất cả các đỉnh trung gian `wᵢ` là mới.

1. **Nền tảng:** Bắt đầu với tam giác `{0, 1, 2}` — đồ thị 2-liên thông đỉnh nhỏ nhất.
2. **Gắn tai:** Trộn ngẫu nhiên các đỉnh còn lại. Lặp lại: lấy 1–3 đỉnh mới, chọn ngẫu nhiên 2 đỉnh đã có làm hai đầu tai, xây dựng đường đi mở. Mỗi tai duy trì tính 2-liên thông đỉnh (định lý Whitney).
3. **Cạnh thêm:** Điền các vị trí `E - |cạnh hiện tại|` còn lại bằng các cặp đỉnh chưa có cạnh.
4. **Hoán vị nhãn:** Áp dụng hoán vị đỉnh ngẫu nhiên.

**Tài liệu tham khảo:** Whitney, H. (1932). *Non-separable and planar graphs.* Trans. AMS, 34(2), 339–362.  
**Độ phức tạp:** O(V²)

---

## Cây vô hướng ngẫu nhiên
**File:** `pycbggp/constructivemethods/undirected_tree.py`  
**Hàm:** `generate_undirected_tree(nb_nodes)`

### Thuật toán — Giải mã dãy Prüfer (Prüfer, 1918)
Dãy Prüfer độ dài `V-2` trên tập `{0..V-1}` mã hóa song ánh một cây có nhãn trên V đỉnh (theo công thức Cayley: có `V^(V-2)` cây có nhãn phân biệt). Sinh dãy ngẫu nhiên đều → thu được **cây có nhãn ngẫu nhiên đều**.

**Quy trình giải mã:**
1. Tính `degree[v]` = số lần xuất hiện của `v` trong dãy Prüfer + 1.
2. Khởi tạo min-heap với tất cả lá (degree == 1).
3. Với mỗi phần tử `p` trong dãy:
   - Lấy lá nhỏ nhất `l` từ heap.
   - Thêm cạnh `(l, p)`, giảm `degree[l]` và `degree[p]`.
   - Nếu `degree[p]` thành 1, đẩy `p` vào heap.
4. Hai đỉnh cuối cùng có degree bằng 1 tạo cạnh cuối.
5. Áp dụng hoán vị đỉnh ngẫu nhiên.

**Tài liệu tham khảo:** Prüfer, H. (1918). *Neuer Beweis eines Satzes über Permutationen.* Arch. Math. Phys., 27, 142–144.  
**Độ phức tạp:** O(V log V)  
**Phân phối:** Đều trên toàn bộ V^(V-2) cây có nhãn.

---

## Cây vô hướng có giới hạn bậc và đường kính
**File:** `pycbggp/constructivemethods/undirected_tree_bounded_diameter_degree.py`  
**Hàm:** `generate_undirected_tree_bounded_diameter_degree(nb_nodes, ub_deg, ub_diameter)`

### Điều kiện khả thi
Đặt `short = ub_diameter // 2`, `long = ub_diameter - short`.  
Kích thước cây tối đa (gốc tại tâm):
- **Đường kính chẵn:** `max_n = 1 + ub_deg × S(short-1, ub_deg-1)`
- **Đường kính lẻ:** `max_n = 1 + S(long-1, ub_deg-1) + (ub_deg-1) × S(short-1, ub_deg-1)`

trong đó `S(d, c) = (c^(d+1) - 1) / (c - 1)` là kích thước cây con tối đa với chiều cao `d` và hệ số nhánh `c`.

### Thuật toán — Gắn đỉnh theo slot có giới hạn độ sâu

**Ý tưởng chính:** Đặt gốc cây tại tâm. Đường đi giữa hai lá đều đi qua gốc, nên `diameter = depth(lá₁) + depth(lá₂)`. Bằng cách giới hạn một nhánh là `long_depth` và tất cả nhánh còn lại là `short_depth`, ta đảm bảo `diameter ≤ long + short = ub_diameter`.

**Mô hình slot:** Duy trì danh sách `avail` gồm các phần tử `[parent, child_depth_budget, slots]`:
- Gốc có 1 slot với ngân sách `long_depth` và `D-1` slot với ngân sách `short_depth`.
- Khi thêm đỉnh `v` vào slot có ngân sách `b`: `v` đóng góp `D-1` slot với ngân sách `b-1` vào `avail` (nếu `b-1 > 0`).

**Mỗi vòng lặp:**
1. Chọn ngẫu nhiên một slot từ `avail` (có trọng số theo số lượng slot).
2. Tính `min_take` và `max_take` dựa trên số đỉnh còn lại và khả năng chứa của các slot.
3. Lấy số lượng con ngẫu nhiên từ slot, gắn vào cây, thêm các slot con của chúng.
4. Giảm số slot; xóa slot trống.
5. Trộn nhãn đỉnh ở cuối.

**Đảm bảo đường kính:** Ràng buộc cứng — không cần retry. Lệnh `assert` ở cuối chỉ là bẫy phát hiện lỗi, không kỳ vọng được kích hoạt.  
**Độ phức tạp:** O(V × D) trong trường hợp xấu nhất.

---

## Đồ thị vô hướng có ràng buộc (V, E, C, B, A)
**File:** `pycbggp/constructivemethods/constrained_graph.py`  
**Các hàm:** `generate_constructive`, `tarjan_analysis`, `check_feasibility`, `verify_graph`

Sinh đồ thị vô hướng thỏa mãn **5 ràng buộc cấu trúc đồng thời**:
- V = số đỉnh, E = số cạnh, C = số thành phần liên thông
- B = số cầu, A = số điểm khớp

### Phân tích Tarjan-Hopcroft — O(V+E)
DFS lặp tính: thành phần liên thông, cầu (qua `low[u] > dfn[parent]`), điểm khớp, và thành phần 2-liên thông cạnh (trích xuất BCC qua stack). Dùng để xác minh kết quả.

### Kiểm tra tính khả thi
Điều kiện cần rút ra từ lý thuyết đồ thị:
- `V-C ≤ E ≤ V_main*(V_main-1)/2` (giới hạn số cạnh)
- `B ≤ V-C` (giới hạn số cầu)
- `A ≤ V_main-2` (giới hạn số điểm khớp)
- Nếu `B=0, A≥1`: cần `V_main ≥ 2A+3` và `E ≥ A+V_main`
- `E-B ∉ {1,2}` (số cạnh không phải cầu phải tạo được chu trình, cần ≥ 3 cạnh)

### Thuật toán Constructive 4 pha
**Pha 1:** `C-1` đỉnh cô lập tạo thêm các thành phần liên thông.

**Pha 2a — Xương sống cầu:**
- `B=1`: một cạnh cầu duy nhất.
- `B≥2, A≤B-1`: đồ thị sâu bướm (sao hoặc đường đi + nhánh treo) tạo đúng A điểm khớp.
- `B≥2, A>B-1`: đường đi thuần B cầu tạo B-1 điểm khớp (bổ sung ở Pha 2b).

**Pha 2b — Khối chu trình (thêm điểm khớp):**
- `B=0, A=0`: chu trình Hamilton (nền tảng 2-liên thông đỉnh).
- `B=0, A≥1`: chuỗi A+1 tam giác.
- `A_extra > 0`: thêm chuỗi tam giác từ các lá của cây cầu (mỗi tam giác thêm 1 điểm khớp).

**Pha 3 — Đệm đỉnh:**
Chia cạnh trong các khối 2-liên thông hiện có: thay cạnh `(u,v)` bằng đường đi `u→x→v`. Thêm đỉnh mà không thay đổi số cầu và điểm khớp.

**Pha 4 — Dày đặc hóa cạnh:**
Điền `E - |hiện tại|` cạnh còn lại trong các khối hiện có (điền hoàn chỉnh đồ thị trong từng khối).

**Độ phức tạp:** O(V+E)  
**Tính đúng đắn:** Chính xác — sinh ra đúng bộ (V,E,C,B,A) cho mọi đầu vào khả thi.

---

## Phương pháp heuristic cho đồ thị có ràng buộc (MCMC, Simulated Annealing)
**File:** `heuristicmethods/constrained_graph.py`  
**Các hàm:** `generate_mcmc`, `generate_simulated_annealing`, `generate_graph`

Bổ sung cho thuật toán constructive: dùng để sinh đồ thị **đa dạng ngẫu nhiên** (MCMC) hoặc làm **phương án dự phòng** cho các trường hợp đặc biệt (SA).

### MCMC — Đổi chỗ cạnh ngẫu nhiên
Bắt đầu từ đồ thị constructive, sau đó áp dụng các phép hoán đổi cạnh đôi ngẫu nhiên trong một thành phần 2-liên thông cạnh. Phép hoán đổi `(a-b, c-d) → (a-c, b-d)` chỉ được chấp nhận nếu (C,B,A) được bảo toàn (kiểm tra qua Tarjan). Sau `nb_iterations` lần hoán đổi, đồ thị giữ nguyên cấu trúc nhưng có topo cạnh khác.

### Simulated Annealing — Luyện thép giả lập
Hàm phạt: `F = w_C*|C_curr-C| + w_B*|B_curr-B| + w_A*|A_curr-A|`.  
Đột biến: xóa 1 cạnh ngẫu nhiên, thêm 1 cạnh ngẫu nhiên mới.  
Chấp nhận Metropolis: trạng thái xấu hơn được chấp nhận với xác suất `exp(-ΔF/T)`, nhiệt độ giảm theo hệ số α mỗi bước.
