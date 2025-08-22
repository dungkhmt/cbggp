#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long
#define dl double
#define st first
#define nd second
#define II pair <int, int>

using namespace std;

const int inf = 7 + 1e9;

class DynamicBitset {
private:
    uint64_t* data;
    int size_;
    int words_;
    static constexpr int BITS_PER_WORD = 64;

    static constexpr uint64_t ALL_ONES = ~0ULL;

    inline int word_index(int pos) const noexcept { return pos >> 6; }
    inline int bit_index(int pos) const noexcept { return pos & 63; }
    
    // Helper function for aligned allocation on Windows
    static uint64_t* allocate_aligned(size_t size) {
#ifdef _WIN32
        return static_cast<uint64_t*>(_aligned_malloc(size, 64));
#else
        return static_cast<uint64_t*>(aligned_alloc(64, size));
#endif
    }
    
    static void deallocate_aligned(uint64_t* ptr) {
#ifdef _WIN32
        _aligned_free(ptr);
#else
        free(ptr);
#endif
    }
    
public:
    // Constructor with memory alignment for better cache performance
    explicit DynamicBitset(int n) noexcept : size_(n) {
        words_ = (n + BITS_PER_WORD - 1) >> 6;
        data = allocate_aligned(words_ * sizeof(uint64_t));
        memset(data, 0, words_ * sizeof(uint64_t));
    }

    DynamicBitset() : size_(0), words_(0), data(nullptr) {}

    DynamicBitset(const DynamicBitset& other) noexcept 
        : size_(other.size_), words_(other.words_) {
        data = allocate_aligned(words_ * sizeof(uint64_t));
        memcpy(data, other.data, words_ * sizeof(uint64_t));
    }

    DynamicBitset(DynamicBitset&& other) noexcept 
        : data(other.data), size_(other.size_), words_(other.words_) {
        other.data = nullptr;
        other.size_ = 0;
        other.words_ = 0;
    }

    DynamicBitset& operator=(const DynamicBitset& other) noexcept {
        if (this != &other) {
            if (words_ != other.words_) {
                deallocate_aligned(data);
                words_ = other.words_;
                data = allocate_aligned(words_ * sizeof(uint64_t));
            }
            size_ = other.size_;
            memcpy(data, other.data, words_ * sizeof(uint64_t));
        }
        return *this;
    }
    
    DynamicBitset& operator=(DynamicBitset&& other) noexcept {
        if (this != &other) {
            deallocate_aligned(data);
            data = other.data;
            size_ = other.size_;
            words_ = other.words_;
            other.data = nullptr;
            other.size_ = 0;
            other.words_ = 0;
        }
        return *this;
    }
    
    ~DynamicBitset() { 
        deallocate_aligned(data); 
    }

    inline void set(int pos) noexcept {
        data[word_index(pos)] |= (1ULL << bit_index(pos));
    }
    
    inline void reset(int pos) noexcept {
        data[word_index(pos)] &= ~(1ULL << bit_index(pos));
    }
    
    inline void flip(int pos) noexcept {
        data[word_index(pos)] ^= (1ULL << bit_index(pos));
    }
    
    inline bool test(int pos) const noexcept {
        return (data[word_index(pos)] >> bit_index(pos)) & 1;
    }

    inline void set_all() noexcept {
        for (int i = 0; i < words_; ++i) {
            data[i] = ALL_ONES;
        }
        int extra = size_ % BITS_PER_WORD;
        if (extra != 0) {
            data[words_ - 1] &= (1ULL << extra) - 1;
        }
    }

    DynamicBitset operator&(const DynamicBitset& other) const noexcept {
        DynamicBitset result(size_);
        const int unroll_factor = 4;
        int i = 0;
    
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            result.data[i] = data[i] & other.data[i];
            result.data[i+1] = data[i+1] & other.data[i+1];
            result.data[i+2] = data[i+2] & other.data[i+2];
            result.data[i+3] = data[i+3] & other.data[i+3];
        }
    
        for (; i < words_; ++i) {
            result.data[i] = data[i] & other.data[i];
        }
        return result;
    }
    
    DynamicBitset operator|(const DynamicBitset& other) const noexcept {
        DynamicBitset result(size_);
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            result.data[i] = data[i] | other.data[i];
            result.data[i+1] = data[i+1] | other.data[i+1];
            result.data[i+2] = data[i+2] | other.data[i+2];
            result.data[i+3] = data[i+3] | other.data[i+3];
        }
        
        for (; i < words_; ++i) {
            result.data[i] = data[i] | other.data[i];
        }
        return result;
    }
    
    DynamicBitset operator^(const DynamicBitset& other) const noexcept {
        DynamicBitset result(size_);
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            result.data[i] = data[i] ^ other.data[i];
            result.data[i+1] = data[i+1] ^ other.data[i+1];
            result.data[i+2] = data[i+2] ^ other.data[i+2];
            result.data[i+3] = data[i+3] ^ other.data[i+3];
        }
        
        for (; i < words_; ++i) {
            result.data[i] = data[i] ^ other.data[i];
        }
        return result;
    }

    DynamicBitset& operator&=(const DynamicBitset& other) noexcept {
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            data[i] &= other.data[i];
            data[i+1] &= other.data[i+1];
            data[i+2] &= other.data[i+2];
            data[i+3] &= other.data[i+3];
        }
        
        for (; i < words_; ++i) {
            data[i] &= other.data[i];
        }
        return *this;
    }
    
    DynamicBitset& operator|=(const DynamicBitset& other) noexcept {
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            data[i] |= other.data[i];
            data[i+1] |= other.data[i+1];
            data[i+2] |= other.data[i+2];
            data[i+3] |= other.data[i+3];
        }
        
        for (; i < words_; ++i) {
            data[i] |= other.data[i];
        }
        return *this;
    }
    
    DynamicBitset& operator^=(const DynamicBitset& other) noexcept {
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            data[i] ^= other.data[i];
            data[i+1] ^= other.data[i+1];
            data[i+2] ^= other.data[i+2];
            data[i+3] ^= other.data[i+3];
        }
        
        for (; i < words_; ++i) {
            data[i] ^= other.data[i];
        }
        return *this;
    }

    inline int count() const noexcept {
        int cnt = 0;
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            cnt += __builtin_popcountll(data[i]);
            cnt += __builtin_popcountll(data[i+1]);
            cnt += __builtin_popcountll(data[i+2]);
            cnt += __builtin_popcountll(data[i+3]);
        }
        
        for (; i < words_; ++i) {
            cnt += __builtin_popcountll(data[i]);
        }
        return cnt;
    }

    inline bool any() const noexcept {
        const int unroll_factor = 4;
        int i = 0;
        
        for (; i <= words_ - unroll_factor; i += unroll_factor) {
            if (data[i] | data[i+1] | data[i+2] | data[i+3]) return true;
        }
        
        for (; i < words_; ++i) {
            if (data[i]) return true;
        }
        return false;
    }
    
    inline bool none() const noexcept { return !any(); }

    template<typename Func>
    inline void for_each_set_bit(Func&& func) const noexcept {
        for (int word_idx = 0; word_idx < words_; ++word_idx) {
            uint64_t word = data[word_idx];
            while (word) {
                int bit = __builtin_ctzll(word);
                func(word_idx * BITS_PER_WORD + bit);
                word &= word - 1;
            }
        }
    }

    vector<int> get_set_bits() const {
        vector<int> result;
        result.reserve(count());
        
        for (int word_idx = 0; word_idx < words_; ++word_idx) {
            uint64_t word = data[word_idx];
            while (word) {
                int bit = __builtin_ctzll(word);
                result.push_back(word_idx * BITS_PER_WORD + bit);
                word &= word - 1;
            }
        }
        return result;
    }

    inline int find_first() const noexcept {
        for (int i = 0; i < words_; ++i) {
            if (data[i]) {
                return i * BITS_PER_WORD + __builtin_ctzll(data[i]);
            }
        }
        return size_;
    }
    
    inline int find_next(int pos) const noexcept {
        if (++pos >= size_) return size_;
        
        int word_idx = word_index(pos);
        int bit_idx = bit_index(pos);

        uint64_t word = data[word_idx] & (~0ULL << bit_idx);
        if (word) {
            return word_idx * BITS_PER_WORD + __builtin_ctzll(word);
        }
    
        for (int i = word_idx + 1; i < words_; ++i) {
            if (data[i]) {
                return i * BITS_PER_WORD + __builtin_ctzll(data[i]);
            }
        }
        return size_;
    }
    
    inline int size() const noexcept { return size_; }
    inline void clear() noexcept { memset(data, 0, words_ * sizeof(uint64_t)); }
    
    inline bool operator[](int pos) const noexcept { return test(pos); }
    
    bool operator==(const DynamicBitset& other) const noexcept {
        if (size_ != other.size_) return false;
        return memcmp(data, other.data, words_ * sizeof(uint64_t)) == 0;
    }
    
    bool operator!=(const DynamicBitset& other) const noexcept {
        return !(*this == other);
    }
};

random_device rd;
mt19937 gen(rd());

int randomInt(const int &l, const int &r) {
    uniform_int_distribution<int> dis(l, r);
    return dis(gen);
}

struct VertexForGreedyAlgorithm {
    int id, degree, frequency, inSolution, tabuTenure;
    VertexForGreedyAlgorithm(const int &_id) {
        id = _id;
        degree = 0;
        frequency = 0;
        inSolution = 0;
        tabuTenure = 0;
    }
};

class GraphForGreedyAlgorithm {
private:
    int numberVertices, numberEdges;
    vector <VertexForGreedyAlgorithm> *vertices;
    vector <II> *edges;
    vector <vector <int>> adjacency, connected;
    dl density;

    vector <int> amts(const int &k, const int &L, const int &iterMax) {
        auto currentSolution = initialize(k);
        updateSolution(currentSolution);
        int iter = 0;
        while (iter < iterMax) {
            auto newSolution = tabuSearch(currentSolution, k, L, iter);
            updateSolution(newSolution);
            int costOfNewSolution = getCostOfSolution(newSolution);
            if (costOfNewSolution * 2 == k * (k - 1))
                return newSolution;
            frequencyBasedInitialize(k, currentSolution);
            updateSolution(currentSolution);
        }
        return vector <int> ();
    }

    vector <int> initialize(const int &k) {
        auto &_vertices = *vertices;

        for (int u = 0; u < numberVertices; u ++) {
            _vertices[u].degree = 0;
            _vertices[u].frequency = 0;
            _vertices[u].inSolution = 0;
            _vertices[u].tabuTenure = 0;
        }
        vector <int> solution;
        for (int i = 0; i < k; i ++) {
            int maxDegree = 0;
            for (int v = 0; v < numberVertices; v ++)
                if (!_vertices[v].inSolution && _vertices[v].degree > maxDegree)
                    maxDegree = _vertices[v].degree;
            vector <int> listOfVertices;
            for (int v = 0; v < numberVertices; v ++)
                if (!_vertices[v].inSolution && _vertices[v].degree == maxDegree)
                    listOfVertices.push_back(v);
            int u = listOfVertices[randomInt(0, (int)listOfVertices.size() - 1)];
            solution.push_back(u);
            _vertices[u].inSolution = 1;
            for (int v : adjacency[u])
                _vertices[v].degree ++;
        }
        return solution;
    }

    vector <int> tabuSearch(vector <int> &solution, const int &k, const int &L, int &iter) {
        int costOfSolution = getCostOfSolution(solution);
        auto bestSolution = solution;
        int costOfBestSolution = costOfSolution, i = 0;
        while (i < L) {
            iter ++;
            auto [u, v, delta] = chooseSwapMove(k, iter, costOfSolution);
            if (u < 0 || v < 0)
                break;
            updateSolutionBySwap(solution, u, v);
            costOfSolution += delta;
            updateTabuList(u, v, k, iter, costOfSolution);
            if (costOfSolution * 2 == k * (k - 1))
                return solution;
            if (costOfSolution > costOfBestSolution) {
                bestSolution = solution;
                costOfBestSolution = costOfSolution;
                i = 0;
            }
            else
                i ++;
        }
        return bestSolution;
    }

    void updateTabuList(const int &u, const int &v, const int &k, const int &iter, const int &costOfSolution) {
        auto &vertices_ = *vertices;

        int l = min(k * (k - 1) / 2 - costOfSolution, 10),
            C = max(k / 40, 6);
        vertices_[u].tabuTenure = iter + l + randomInt(0, C);
        vertices_[v].tabuTenure = iter + int(l * 0.6) + randomInt(0, int(C * 0.6));
    }

    void updateSolutionBySwap(vector <int> &solution, const int &u, const int &v) {
        auto &_vertices = *vertices;

        solution.erase(remove(solution.begin(), solution.end(), u), solution.end());
        _vertices[u].inSolution = 0;
        for (int w : adjacency[u])
            _vertices[w].degree --;
        _vertices[u].frequency ++;

        solution.push_back(v);
        _vertices[v].inSolution = 1;
        for (int w : adjacency[v])
            _vertices[w].degree ++;
        _vertices[v].frequency ++;
    }

    tuple <int, int, int> chooseSwapMove(const int &k, const int &iter, const int &costOfSolution) {
        auto &vertices_ = *vertices;

        vector <int> A, B;
        int minInSolution = k, maxOutSolution = 0;

        for (int u = 0; u < numberVertices; u ++)
            if (vertices_[u].tabuTenure < iter) {
                if (vertices_[u].inSolution)
                    minInSolution = min(minInSolution, vertices_[u].degree);
                else
                    maxOutSolution = max(maxOutSolution, vertices_[u].degree);
            }

        for (int u = 0; u < numberVertices; u ++)
            if (vertices_[u].tabuTenure < iter) {
                if (vertices_[u].inSolution && vertices_[u].degree == minInSolution)
                    A.push_back(u);
                else if (!vertices_[u].inSolution && vertices_[u].degree == maxOutSolution)
                    B.push_back(u);
            }

        if (A.empty() || B.empty())
            return {-1, -1, 0};

        auto T = bestSwapMoves(A, B);
        int l = k * (k - 1) / 2 - costOfSolution;
        dl p = min((l + 2.0) / numberEdges, 0.1),
           value = static_cast<double>(rand()) / RAND_MAX;

        if ((maxOutSolution <= minInSolution || (T.empty() && maxOutSolution == minInSolution + 1)) && value <= p)
            return pmsr(k, iter);

        if (!T.empty()) {
            auto [u, v] = T[randomInt(0, (int)T.size() - 1)];
            return {u, v, maxOutSolution - minInSolution};
        }
        else {
            int u = A[randomInt(0, (int)A.size() - 1)],
                v = B[randomInt(0, (int)B.size() - 1)];
            return {u, v, maxOutSolution - minInSolution - 1};
        }
    }
    
    vector <II> bestSwapMoves(const vector <int> &A, const vector <int> &B) {
        if (A.empty() || B.empty())
            return vector <II> ();
        vector <II> T;
        for (int u : A)
            for (int v : B)
                if (!connected[u][v])
                    T.push_back({u, v});
        return T;
    }

    tuple <int, int, int> pmsr(const int &k, const int &iter) {
        auto &vertices_ = *vertices;

        vector <int> listOfVertices;
        for (int v = 0; v < numberVertices; v ++)
            listOfVertices.push_back(v);
        shuffle(listOfVertices.begin(), listOfVertices.end(), gen);

        int u = -1, v = -1;
        for (int w : listOfVertices) {
            if (vertices_[w].tabuTenure >= iter)
                continue;
            if (u < 0 && vertices_[w].inSolution)
                u = w;
            if (v < 0 && !vertices_[w].inSolution && vertices_[w].degree < int(k * density))
                v = w;
            if (u >= 0 && v >= 0)
                break;
        }
        return {u, v, vertices_[v].degree - vertices_[u].degree - (connected[u][v] ? 1 : 0)};
    }

    void updateSolution(vector <int> &solution) {
        auto &vertices_ = *vertices;

        for (int id = 0; id < numberVertices; id ++) {
            vertices_[id].degree = 0;
            vertices_[id].frequency = 0;
            vertices_[id].inSolution = 0;
            vertices_[id].tabuTenure = 0;
        }

        for (int u : solution) {
            vertices_[u].inSolution = 1;
            for (int v : adjacency[u])
                vertices_[v].degree ++;
        }
    }

    int getCostOfSolution(const vector <int> &solution) {
        int cost = 0;
        for (int u : solution)
            cost += vertices->at(u).degree;
        return cost / 2;
    }

    vector <int> frequencyBasedInitialize(const int &k, vector <int> &solution) {
        auto &vertices_ = *vertices;

        refreshSolution(k);
        solution.clear();
        int minFrequency = numberVertices;
        for (const auto &v : vertices_)
            if (v.frequency < minFrequency)
                minFrequency = v.frequency;

        vector <int> listVerticesMinFrequency;
        for (int i = 0; i < numberVertices; i ++)
            if (vertices_[i].frequency == minFrequency)
                listVerticesMinFrequency.push_back(i);

        int u = listVerticesMinFrequency[randomInt(0, (int)listVerticesMinFrequency.size() - 1)];
        solution.push_back(u);
        vertices_[u].inSolution = 1;
        for (int v : adjacency[u])
            vertices_[v].degree ++;

        while ((int)solution.size() < k) {
            int minDegree = numberVertices, minFrequencyOfMinDegree = numberVertices;
            for (int v = 0; v < numberVertices; v ++)
                if (!vertices_[v].inSolution) {
                    if (vertices_[v].degree < minDegree) {
                        minDegree = vertices_[v].degree;
                        minFrequencyOfMinDegree = vertices_[v].frequency;
                    }
                    else if (vertices_[v].degree == minDegree && vertices_[v].frequency < minFrequencyOfMinDegree)
                        minFrequencyOfMinDegree = vertices_[v].frequency;
                }

            vector <int> listOfVertices;
            for (int v = 0; v < numberVertices; v ++)
                if (!vertices_[v].inSolution && vertices_[v].degree == minDegree && vertices_[v].frequency == minFrequencyOfMinDegree)
                    listOfVertices.push_back(v);

            u = listOfVertices[randomInt(0, (int)listOfVertices.size() - 1)];
            solution.push_back(u);
            vertices_[u].inSolution = 1;
            for (int v : adjacency[u])
                vertices_[v].degree ++;
        }
        return solution;
    }

    void refreshSolution(const int &k) {
        auto &vertices_ = *vertices;

        for (int id = 0; id < numberVertices; id ++) {
            vertices_[id].degree = 0;
            vertices_[id].frequency %= k;
            vertices_[id].inSolution = 0;
        }
    }

public:
    GraphForGreedyAlgorithm(const int &_numberVertices, vector <VertexForGreedyAlgorithm> &_vertices, vector <II> &_edges) {
        numberVertices = _numberVertices;
        numberEdges = _edges.size();
        vertices = &_vertices;
        edges = &_edges;
        adjacency.resize(numberVertices);
        connected.assign(numberVertices, vector<int>(numberVertices, 0));
        for (const auto& edge : *edges) {
            int u = edge.st, v = edge.nd;
            adjacency[u].push_back(v);
            adjacency[v].push_back(u);
            connected[u][v] = 1;
            connected[v][u] = 1;
        }
        density = (numberVertices > 1 ? (dl)numberEdges / (numberVertices * (numberVertices - 1) / 2) : 0);
    }

    vector <int> findClique() {
        vector <int> solution;
        int k = 1;
        while (k <= numberVertices) {
            vector <int> newSolution = amts(k, k * 4, k * 10);
            if ((int)newSolution.size() == k) {
                solution = newSolution;
                k ++;
            }
            else
                break;
        }
        return solution;
    }
};

class Vertex {
public:
    int id, positionInListSorted;

    Vertex() : id(-1), positionInListSorted(-1) {}

    Vertex (const int &_id) : id(_id), positionInListSorted(_id) {}
};

class Graph {
public:
    int numberVertices;
    vector <Vertex> *vertices;
    vector <DynamicBitset> *adjacency;
    vector <vector <int>> connected;

    Graph(const int &_numberVertices, vector <Vertex> &_vertices, vector <DynamicBitset> &_adjacency) {
        numberVertices = _numberVertices;
        vertices = &_vertices;
        adjacency = &_adjacency;
    }

    int checkSolution(vector <int> &solution) {
        if (solution.empty())
            return 0;
        for (int i = 0; i < solution.size(); i ++)
            for (int j = i + 1; j < solution.size(); j ++)
                if (!(*adjacency)[solution[i]].test(solution[j]))
                    return 0;
        return 1;
    }

    vector <int> listVerticesOfNeighbors(const DynamicBitset &adjacencyOfU) {
        vector <int> listOfVertices;
        for (int v = adjacencyOfU.find_first(); v < adjacencyOfU.size(); v = adjacencyOfU.find_next(v))
            listOfVertices.push_back(v);
        return listOfVertices;
    }
};

void degSort(Graph &graph, vector <II> &edges) {
    vector <int> degree(graph.numberVertices, 0);
    for (const auto &edge : edges) {
        int u = edge.st, v = edge.nd;
        degree[u] ++;
        degree[v] ++;
    }
    sort(graph.vertices->begin(), graph.vertices->end(), [&](const Vertex &a, const Vertex &b) {
        if (degree[a.id] != degree[b.id])
            return degree[a.id] > degree[b.id];
        return a.id < b.id;
    });
    auto &vertices = *graph.vertices;
    auto &adjacency = *graph.adjacency;
    vector <int> positions(graph.numberVertices, 0);
    for (int i = 0; i < graph.numberVertices; i ++) {
        vertices[i].positionInListSorted = i;
        positions[vertices[i].id] = i;
    }
    adjacency.assign(graph.numberVertices, DynamicBitset(graph.numberVertices));
    graph.connected.assign(graph.numberVertices, vector <int> (graph.numberVertices, 0));
    for (int i = 0; i < edges.size(); i ++) {
        int u = positions[edges[i].st], v = positions[edges[i].nd];
        edges[i] = {u, v};
        adjacency[u].set(v);
        adjacency[v].set(u);
        graph.connected[u][v] = 1;
        graph.connected[v][u] = 1;
    }
}

vector <int> buildInitialUpperBound(Graph &graph) {
    vector <int> upperBound(graph.numberVertices, 1);
    for (int u = 0; u < graph.numberVertices; u ++) {
        DynamicBitset &adjacencyOfU = graph.adjacency->at(u);
        for (int v = adjacencyOfU.find_first(); v < adjacencyOfU.size(); v = adjacencyOfU.find_next(v)) {
            if (v >= u)
                break;
            upperBound[u] = max(upperBound[u], upperBound[v] + 1);
        }
    }
    return upperBound;
}

class GraphBitString {
private:
    vector <vector <int>> *connected;

public:
    int numberVertices;
    vector <int> vertices;
    vector <DynamicBitset> adjacency;
    unordered_map <int, int> positionInVertices;

    GraphBitString(const int &_numberVertices, vector <vector <int>> &_connected, vector<int> _vertices) {
        numberVertices = _numberVertices;
        connected = &_connected;
        vertices = _vertices;
        sort(vertices.begin(), vertices.end());
        for (int indexOfU = 0; indexOfU < (int)vertices.size(); indexOfU ++)
            positionInVertices[vertices[indexOfU]] = indexOfU;
        adjacency.assign(vertices.size(), DynamicBitset(vertices.size()));
        for (int indexOfU = 0; indexOfU < (int)vertices.size(); indexOfU ++) {
            int u = vertices[indexOfU];
            for (int indexOfV = indexOfU + 1; indexOfV < (int)vertices.size(); indexOfV ++) {
                int v = vertices[indexOfV];
                if (!_connected[u][v])
                    adjacency[indexOfU].set(indexOfV);
            }
        }
    }

    tuple <int, vector <int>> iseq(const int &numberColors, vector <int> &listOfVertices) {
        vector <int> colorOfVertices(numberVertices, -1);

        DynamicBitset bitStringOfVertices(vertices.size());
        for (int u : listOfVertices)
            bitStringOfVertices.set(positionInVertices[u]);
        
        for (int color = 0; color < numberColors; color ++) {
            if (!bitStringOfVertices.any())
                return {0, vector <int> ()};
            
            DynamicBitset copyBitStringOfVertices = bitStringOfVertices;
            while (copyBitStringOfVertices.any()) {
                int indexOfU = copyBitStringOfVertices.find_first();
                colorOfVertices[vertices[indexOfU]] = color;
                bitStringOfVertices.reset(indexOfU);
                copyBitStringOfVertices &= adjacency[indexOfU];
            }
        }
        if (!bitStringOfVertices.any())
            return {0, vector <int> ()};
        return {1, colorOfVertices};
    }

    int filt(const int &numberColors, const int &vertex, vector <int> &lastVertexOfColorHasVertex, vector <int> &verticesOfNewSubGraph, vector <int> &isInNewPrunedSet) {
        DynamicBitset bitStringOfVertices(vertices.size());
        for (int u : verticesOfNewSubGraph)
            if (isInNewPrunedSet[u] && (*connected)[vertex][u])
                bitStringOfVertices.set(positionInVertices[u]);
        
        for (int color = 0; color < numberColors; color ++) {
            if (!bitStringOfVertices.any())
                return 0;

            DynamicBitset copyBitStringOfVertices = bitStringOfVertices;
            int indexOfU = copyBitStringOfVertices.find_first(),
                lastVertex = lastVertexOfColorHasVertex[vertices[indexOfU]];
            bitStringOfVertices.reset(indexOfU);
            copyBitStringOfVertices &= adjacency[indexOfU];

            while (copyBitStringOfVertices.any()) {
                indexOfU = copyBitStringOfVertices.find_first();
                int u = vertices[indexOfU];

                if (u < lastVertex)
                    copyBitStringOfVertices &= adjacency[indexOfU];
                else {
                    isInNewPrunedSet[u] = 0;
                    copyBitStringOfVertices.reset(indexOfU);
                }

                bitStringOfVertices.reset(indexOfU);
            }
        }
        return 1;
    }
};

struct SoftClause {
    DynamicBitset bitStringOfVertices;
    int state, numberOfUndefinedLiterals, isInProperSet;

    SoftClause(const DynamicBitset &_bitStringOfVertices) {
        bitStringOfVertices = _bitStringOfVertices;
        state = 0;
        numberOfUndefinedLiterals = 0;
        isInProperSet = 0;
    }
};

class GraphSAT {
private:
    int maxNumber, numberVertices;
    vector <int> vertices, stateOfVertex;
    vector <DynamicBitset> adjacency;
    
    void updateStateOfVertex(const int &indexOfVertex, const int &state, const int &indexOfSoftClause) {
        stateOfVertex[indexOfVertex] = state;

        if (indexOfSoftClause >= 0) {
            softClauses[indexOfSoftClause].state += state;
            softClauses[indexOfSoftClause].numberOfUndefinedLiterals -= 1;
            if (state == 1)
                softClauses[indexOfSoftClause].isInProperSet = 1;
        }
    }

    void addNewVertexToSoftClause(const int &softClauseIndex) {
        int u = maxNumber;
        positionInVertices[u] = numberVertices;
        adjacency.push_back(DynamicBitset(maxSize));
        indexOfSoftClause[numberVertices] = softClauseIndex;

        SoftClause &softClause = softClauses[softClauseIndex];        
        DynamicBitset &bitStringOfVertices = softClause.bitStringOfVertices;
        for (int indexOfV = bitStringOfVertices.find_first(); indexOfV < bitStringOfVertices.size(); indexOfV = bitStringOfVertices.find_next(indexOfV)) {
            adjacency[indexOfV].set(numberVertices);
            adjacency[numberVertices].set(indexOfV);
        }
        
        softClauses[softClauseIndex].bitStringOfVertices.set(numberVertices);

        numberVertices ++;
        maxNumber ++;
    }

public:
    int maxSize;
    unordered_map <int, int> positionInVertices, indexOfSoftClause;
    vector <SoftClause> softClauses;

    GraphSAT(GraphBitString &bitStringSubGraph, vector <int> _vertices) {
        maxNumber = bitStringSubGraph.numberVertices;
        numberVertices = _vertices.size();
        maxSize = numberVertices * 2;
        vertices = _vertices;
        DynamicBitset bitStringOfVertices(maxNumber);
        for (int u : vertices)
            bitStringOfVertices.set(bitStringSubGraph.positionInVertices[u]);
        for (int indexOfU = 0; indexOfU < numberVertices; indexOfU ++) {
            positionInVertices[vertices[indexOfU]] = indexOfU;
            indexOfSoftClause[indexOfU] = -1;
        }
        adjacency.assign(numberVertices, DynamicBitset(maxSize));
        for (int indexOfU = 0; indexOfU < numberVertices; indexOfU ++) {
            int u = vertices[indexOfU];
            DynamicBitset adjacencyOfU = bitStringSubGraph.adjacency[bitStringSubGraph.positionInVertices[u]] & bitStringOfVertices;
            for (int indexOfVInBitStringSubGraph = adjacencyOfU.find_first(); indexOfVInBitStringSubGraph < adjacencyOfU.size(); indexOfVInBitStringSubGraph = adjacencyOfU.find_next(indexOfVInBitStringSubGraph)) {
                int indexOfV = positionInVertices[bitStringSubGraph.vertices[indexOfVInBitStringSubGraph]];
                adjacency[indexOfU].set(indexOfV);
                adjacency[indexOfV].set(indexOfU);
            }
        }
    }

    void addSoftClause(DynamicBitset &bitStringOfVertices) {
        if (!bitStringOfVertices.any())
            return;

        int clauseIndex = softClauses.size();
        softClauses.push_back(SoftClause(bitStringOfVertices));
        for (int indexOfU = bitStringOfVertices.find_first(); indexOfU < bitStringOfVertices.size(); indexOfU = bitStringOfVertices.find_next(indexOfU))
            indexOfSoftClause[indexOfU] = clauseIndex;
    }

    int unitPropagation(const int &vertex) {
        stateOfVertex.assign(numberVertices, -1);

        for (auto &softClause : softClauses) {
            softClause.state = 0;
            softClause.numberOfUndefinedLiterals = softClause.bitStringOfVertices.count();
            softClause.isInProperSet = 0;
        }

        int indexOfVertex = positionInVertices[vertex];
        stateOfVertex[indexOfVertex] = 1;

        stack <II> s;
        s.push({indexOfVertex, 1});
        softClauses[indexOfSoftClause[indexOfVertex]].isInProperSet = 1;

        while (!s.empty()) {
            auto [indexOfU, state] = s.top();
            s.pop();

            if (state == 1) {
                auto &adjacencyOfU = adjacency[indexOfU];
                for (int indexOfV = adjacencyOfU.find_first(); indexOfV < adjacencyOfU.size(); indexOfV = adjacencyOfU.find_next(indexOfV))
                    if (stateOfVertex[indexOfV] == -1) {
                        updateStateOfVertex(indexOfV, 0, indexOfSoftClause[indexOfV]);
                        s.push({indexOfV, 0});
                    }
            }
            if (indexOfSoftClause[indexOfU] >= 0) {
                auto &softClause = softClauses[indexOfSoftClause[indexOfU]];
                auto &bitStringOfSoftClause = softClause.bitStringOfVertices;

                if (state == 1 && softClause.numberOfUndefinedLiterals > 0)
                    for (int indexOfV = bitStringOfSoftClause.find_first(); indexOfV < bitStringOfSoftClause.size(); indexOfV = bitStringOfSoftClause.find_next(indexOfV))
                        if (stateOfVertex[indexOfV] == -1) {
                            updateStateOfVertex(indexOfV, 0, indexOfSoftClause[indexOfV]);
                            s.push({indexOfV, 0});
                        }

                if (softClause.state == 0) {
                    if (softClause.numberOfUndefinedLiterals == 0) {
                        softClause.isInProperSet = 1;
                        return 0;
                    }
                    if (softClause.numberOfUndefinedLiterals == 1)
                        for (int indexOfV = bitStringOfSoftClause.find_first(); indexOfV < bitStringOfSoftClause.size(); indexOfV = bitStringOfSoftClause.find_next(indexOfV))
                            if (stateOfVertex[indexOfV] == -1) {
                                updateStateOfVertex(indexOfV, 1, indexOfSoftClause[indexOfV]);
                                s.push({indexOfV, 1});
                                break;
                            }
                }
            }
        }
        
        return 1;
    }

    void removeVertex(const int &u) {
        int indexOfU = positionInVertices[u];
        auto &adjacencyOfU = adjacency[indexOfU];
        for (int indexOfV = adjacencyOfU.find_first(); indexOfV < adjacencyOfU.size(); indexOfV = adjacencyOfU.find_next(indexOfV))
            adjacency[indexOfV].reset(indexOfU);
        
        int clauseIndex = indexOfSoftClause[indexOfU];
        if (clauseIndex >= 0) {
            indexOfSoftClause[indexOfU] = -1;
            softClauses[clauseIndex].bitStringOfVertices.reset(indexOfU);
        }
    }

    void removeLastSoftClause() {
        if (softClauses.empty())
            return;
        auto &lastClause = softClauses.back();
        softClauses.pop_back();

        auto &bitStringOfVertices = lastClause.bitStringOfVertices;
        for (int indexOfU = bitStringOfVertices.find_first(); indexOfU < bitStringOfVertices.size(); indexOfU = bitStringOfVertices.find_next(indexOfU))
            indexOfSoftClause[indexOfU] = -1;
    }

    void transformGraph() {
        for (int softClauseIndex = 0; softClauseIndex < (int)softClauses.size(); softClauseIndex ++)
            if (softClauses[softClauseIndex].isInProperSet)
                addNewVertexToSoftClause(softClauseIndex);
    }
};

int filtcol(GraphBitString &bitStringSubGraph, const int &numberColors, vector <int> &lastVertexOfColorHasVertex, vector <int> &verticesOfNewSubGraph, vector <int> &isInNewPrunedSet, vector <int> &isInNewBranchingSet) {
    for (int vertex : verticesOfNewSubGraph) {
        if (!isInNewBranchingSet[vertex])
            continue;
        
        if (!bitStringSubGraph.filt(numberColors - 1, vertex, lastVertexOfColorHasVertex, verticesOfNewSubGraph, isInNewPrunedSet))
            return 0;
    }
    
    return 1;
}

int filtsat(GraphBitString &bitStringSubGraph, const int &numberColors, vector <int> &colorOfVertices, vector <int> &verticesOfNewSubGraph, vector <int> &isInNewPrunedSet, vector <int> &isInNewBranchingSet) {
    GraphSAT subGraphSAT(bitStringSubGraph, verticesOfNewSubGraph);

    vector <DynamicBitset> independentSet(numberColors, DynamicBitset(subGraphSAT.maxSize));
    for (int u : verticesOfNewSubGraph)
        if (isInNewPrunedSet[u])
            independentSet[colorOfVertices[u]].set(subGraphSAT.positionInVertices[u]);
        else if (isInNewBranchingSet[u])
            independentSet[numberColors - 1].set(subGraphSAT.positionInVertices[u]);
    
    for (int color = 0; color < numberColors; color ++)
        subGraphSAT.addSoftClause(independentSet[color]);
    
    for (int u : verticesOfNewSubGraph) {
        if (!isInNewPrunedSet[u] && !isInNewBranchingSet[u])
            continue;
        if (!subGraphSAT.unitPropagation(u)) {
            int indexOfU = subGraphSAT.positionInVertices[u];
            SoftClause &softClause = subGraphSAT.softClauses[subGraphSAT.indexOfSoftClause[indexOfU]];
            if (softClause.bitStringOfVertices.count() == 1)
                return 0;
            
            subGraphSAT.removeVertex(u);
            if (isInNewPrunedSet[u])
                isInNewPrunedSet[u] = 0;
            else
                isInNewBranchingSet[u] = 0;
        }
    }
    
    return 1;
}

void satcol(GraphBitString &bitStringSubGraph, const int &numberColors, vector <int> &colorOfVertices, vector <int> &verticesOfNewSubGraph, vector <int> &isInNewPrunedSet, vector <int> &isInNewBranchingSet) {
    GraphSAT subGraphSAT(bitStringSubGraph, verticesOfNewSubGraph);

    vector <DynamicBitset> independentSet(numberColors, DynamicBitset(subGraphSAT.maxSize));
    DynamicBitset bitStringOfListOfVertices(bitStringSubGraph.vertices.size());
    for (int u : verticesOfNewSubGraph)
        if (isInNewPrunedSet[u])
            independentSet[colorOfVertices[u]].set(subGraphSAT.positionInVertices[u]);
        else
            bitStringOfListOfVertices.set(bitStringSubGraph.positionInVertices[u]);

    for (int color = 0; color < numberColors - 1; color ++)
        subGraphSAT.addSoftClause(independentSet[color]);

    while (bitStringOfListOfVertices.any()) {
        vector <int> independentSetVertices;
        DynamicBitset bitStringOfIndependentSet(subGraphSAT.maxSize);
        DynamicBitset copyBitStringOfListOfVertices = bitStringOfListOfVertices;

        while (copyBitStringOfListOfVertices.any()) {
            int indexOfU = copyBitStringOfListOfVertices.find_first(),
                u = bitStringSubGraph.vertices[indexOfU];
            
            independentSetVertices.push_back(u);
            bitStringOfIndependentSet.set(subGraphSAT.positionInVertices[u]);
            bitStringOfListOfVertices.reset(indexOfU);
            copyBitStringOfListOfVertices &= bitStringSubGraph.adjacency[indexOfU];
        }

        subGraphSAT.addSoftClause(bitStringOfIndependentSet);
        int ok = 0;
        for (int u : independentSetVertices)
            if (subGraphSAT.unitPropagation(u)) {
                ok = 1;
                break;
            }

        if (ok)
            subGraphSAT.removeLastSoftClause();
        else {
            for (int u : independentSetVertices) {
                isInNewPrunedSet[u] = 1;
                isInNewBranchingSet[u] = 0;
            }
            subGraphSAT.transformGraph();
        }
    }
}

int lowerBound = 0;
vector <int> bestSolution, currentSolution;

void findMaxClique(Graph &graph, vector <int> &currentPrunedSet, vector <int> &currentBranchingSet, vector <int> &currentUpperBound) {
    int numberVerticesOfGraph = graph.numberVertices;
    vector <int> newUpperBound = currentUpperBound,
                 isInCurrentPrunedSet(numberVerticesOfGraph, 0);
    
    DynamicBitset bitStringOfCurrentSubGraph(numberVerticesOfGraph);
    for (int u : currentPrunedSet) {
        isInCurrentPrunedSet[u] = 1;
        bitStringOfCurrentSubGraph.set(u);
    }
    for (int u : currentBranchingSet)
        bitStringOfCurrentSubGraph.set(u);

    vector <int> currentSubGraph = currentPrunedSet;
    for (int u : currentBranchingSet)
        currentSubGraph.push_back(u);
    GraphBitString bitStringSubGraph(numberVerticesOfGraph, graph.connected, currentSubGraph);

    auto &adjacency = *graph.adjacency;
    for (int u : currentBranchingSet) {
        vector <int> adjacencyOfU = graph.listVerticesOfNeighbors(adjacency[u] & bitStringOfCurrentSubGraph),
                     verticesOfNewSubGraph;
        for (int v : adjacencyOfU)
            if (isInCurrentPrunedSet[v] || v < u)
                verticesOfNewSubGraph.push_back(v);
        
        newUpperBound[u] = 1;
        for (int v : verticesOfNewSubGraph)
            newUpperBound[u] = max(newUpperBound[u], newUpperBound[v] + 1);
        
        if (newUpperBound[u] + (int)currentSolution.size() <= lowerBound)
            isInCurrentPrunedSet[u] = 1;
        else {
            if (verticesOfNewSubGraph.empty()) {
                if ((int)currentSolution.size() >= (int)bestSolution.size()) {
                    bestSolution = currentSolution;
                    bestSolution.push_back(u);
                    lowerBound = max(lowerBound, (int)bestSolution.size());
                }
                newUpperBound[u] = min(newUpperBound[u], lowerBound - (int)currentSolution.size());
                continue;
            }

            int numberColors = lowerBound - (int)currentSolution.size();
            vector <int> isInNewPrunedSet(numberVerticesOfGraph, 0),
                         isInNewBranchingSet(numberVerticesOfGraph, 1);
            
            if (numberColors > 1) {
                auto [check, colorOfVertices] = bitStringSubGraph.iseq(numberColors - 1, verticesOfNewSubGraph);

                if (!check) {
                    newUpperBound[u] = min(newUpperBound[u], lowerBound - (int)currentSolution.size());
                    continue;
                }

                int isIndependent = 1;
                DynamicBitset bitStringOfBranchingSet(bitStringSubGraph.vertices.size());
                for (int v : verticesOfNewSubGraph)
                    if (colorOfVertices[v] >= 0) {
                        isInNewPrunedSet[v] = 1;
                        isInNewBranchingSet[v] = 0;
                    } else {
                        int indexOfV = bitStringSubGraph.positionInVertices[v];
                        DynamicBitset _ = bitStringSubGraph.adjacency[indexOfV] & bitStringOfBranchingSet;
                        if (_ != bitStringOfBranchingSet)
                            isIndependent = 0;
                        bitStringOfBranchingSet.set(indexOfV);
                    }

                // if (isIndependent) {
                //     vector <int> lastVertexOfColor(numberColors, -1);
                //     for (int v : verticesOfNewSubGraph)
                //         if (isInNewPrunedSet[v])
                //             lastVertexOfColor[colorOfVertices[v]] = v;
                    
                //     vector <int> lastVertexOfColorHasVertex(numberVerticesOfGraph, -1);
                //     for (int v : verticesOfNewSubGraph)
                //         if (isInNewPrunedSet[v])
                //             lastVertexOfColorHasVertex[v] = lastVertexOfColor[colorOfVertices[v]];

                //     if (!filtcol(bitStringSubGraph, numberColors, lastVertexOfColorHasVertex, verticesOfNewSubGraph, isInNewPrunedSet, isInNewBranchingSet) || \
                //        !filtsat(bitStringSubGraph, numberColors, colorOfVertices, verticesOfNewSubGraph, isInNewPrunedSet ,isInNewBranchingSet)) {
                //         newUpperBound[u] = min(newUpperBound[u], lowerBound - (int)currentSolution.size());
                //         continue;
                //     }
                // }
                // else
                    satcol(bitStringSubGraph, numberColors, colorOfVertices, verticesOfNewSubGraph, isInNewPrunedSet, isInNewBranchingSet);
            }

            vector <int> newPrunedSet, newBranchingSet;
            for (int v : verticesOfNewSubGraph) {
                if (isInNewPrunedSet[v])
                    newPrunedSet.push_back(v);
                if (isInNewBranchingSet[v])
                    newBranchingSet.push_back(v);
            }

            if (!newBranchingSet.empty()) {
                currentSolution.push_back(u);
                findMaxClique(graph, newPrunedSet, newBranchingSet, newUpperBound);
                currentSolution.pop_back();
            }
            else if ((int)currentSolution.size() >= (int)bestSolution.size()) {
                bestSolution = currentSolution;
                bestSolution.push_back(u);
                lowerBound = max(lowerBound, (int)bestSolution.size());
            }
            newUpperBound[u] = min(newUpperBound[u], lowerBound - (int)currentSolution.size());
        }
    }
}

int main() {
#define TASKNAME ""
    ios_base :: sync_with_stdio (0);
    cin.tie (0);
    if ( fopen( TASKNAME".inp", "r" ) ) {
        freopen (TASKNAME".inp", "r", stdin);
        freopen (TASKNAME".out", "w", stdout);
    }

    int numberVertices, numberEdges;
    cin >> numberVertices >> numberEdges;
    vector <Vertex> vertices(numberVertices);
    for (int i = 0; i < numberVertices; i ++)
        vertices[i] = Vertex(i);
    vector <II> edges;
    vector <DynamicBitset> adjacency(numberVertices, DynamicBitset(numberVertices));
    for (int i = 0; i < numberEdges; i ++) {
        // char c;
        // cin >> c;
        int u, v;
        cin >> u >> v;
        u --; v --;
        edges.push_back({u, v});
        adjacency[u].set(v);
        adjacency[v].set(u);
    }

    auto startTime = chrono::high_resolution_clock::now();

    Graph graph(numberVertices, vertices, adjacency);
    degSort(graph, edges);

    vector <VertexForGreedyAlgorithm> verticesOfGraphForGreedy;
    for (const auto& vertex : (*graph.vertices))
        verticesOfGraphForGreedy.push_back(VertexForGreedyAlgorithm(vertex.positionInListSorted));
    GraphForGreedyAlgorithm graphForGreedy(numberVertices, verticesOfGraphForGreedy, edges);
    vector <int> initialSolution = graphForGreedy.findClique();
    // vector <int> initialSolution = {0};
    lowerBound = initialSolution.size();
    bestSolution = initialSolution;
    vector <int> initialUpperBound = buildInitialUpperBound(graph);
    
    for (int i = lowerBound; i < numberVertices; i ++) {
        vector <int> verticesOfSubGraph;
        DynamicBitset adjacencyOfI = adjacency[i];
        for (int j = adjacencyOfI.find_first(); j < adjacencyOfI.size(); j = adjacencyOfI.find_next(j)) {
            if (j >= i)
                break;
            verticesOfSubGraph.push_back(j);
        }
        
        if ((int)verticesOfSubGraph.size() < lowerBound)
            continue;

        vector <int> prunedSet(verticesOfSubGraph.begin(), verticesOfSubGraph.begin() + lowerBound),
                     branchingSet(verticesOfSubGraph.begin() + lowerBound, verticesOfSubGraph.end());
        
        currentSolution.push_back(i);
        findMaxClique(graph, prunedSet, branchingSet, initialUpperBound);
        currentSolution.pop_back();
        initialUpperBound[i] = lowerBound;
    }

    auto endTime = chrono::high_resolution_clock::now();
    auto totalTime = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);

    cout << initialSolution.size() << '\n';
    for (int vertex : initialSolution)
        cout << (*graph.vertices)[vertex].id + 1 << ' ';
    cout << "\n\n";
    cout << bestSolution.size() << '\n';
    for (int vertex : bestSolution)
        cout << (*graph.vertices)[vertex].id + 1 << ' ';
    cout << "\n\nCheck Solution: " << (graph.checkSolution(bestSolution) ? "True" : "False") << '\n';
    cout << "Total Time: " << totalTime.count() << " ms";

    return 0;
}
