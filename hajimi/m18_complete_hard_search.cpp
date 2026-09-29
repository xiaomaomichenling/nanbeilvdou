#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif

namespace fs = std::filesystem;
using U = uint32_t;
using Counts = std::map<int, uint64_t>;
constexpr int N = 7;
int M = 18;
U S = U(1) << 18;
constexpr const char *ROOT = "C:/Users/L/Desktop/hajimi";
fs::path OUTPUT_DIR;

struct Pair { int u, v; };
std::array<Pair, 21> pairs;
std::vector<std::array<int, 21>> perm_indices;
std::array<int, 63> cuts;

void init_tables() {
    int p = 0;
    for (int u = 0; u < N; ++u) for (int v = u + 1; v < N; ++v) pairs[p++] = {u, v};
    int c = 0;
    for (int x = 1; x < 127; ++x) if (x & 1) cuts[c++] = x;
    if (c != 63) throw std::runtime_error("cut table failure");
}

using Degrees = std::array<int, N>;
void build_permutations(const Degrees &degrees) {
    perm_indices.clear();
    std::array<int, N> a{0, 1, 2, 3, 4, 5, 6};
    do {
        bool preserves_degrees = true;
        for (int i = 0; i < N; ++i)
            if (degrees[i] != degrees[a[i]]) { preserves_degrees = false; break; }
        if (!preserves_degrees) continue;
        std::array<int, 21> indices{};
        for (int i = 0; i < 21; ++i) {
            int u = a[pairs[i].u], v = a[pairs[i].v];
            if (u > v) std::swap(u, v);
            int j = 0;
            while (pairs[j].u != u || pairs[j].v != v) ++j;
            indices[i] = j;
        }
        perm_indices.push_back(indices);
    } while (std::next_permutation(a.begin(), a.end()));
    if (perm_indices.empty()) throw std::runtime_error("empty degree-preserving permutation group");
}

void degree_sequences_rec(int pos, int ceiling, int sum_left, Degrees &d,
                          std::vector<Degrees> &sequences) {
    if (pos == N) {
        if (sum_left == 0) sequences.push_back(d);
        return;
    }
    int vertices_left = N - pos - 1;
    int upper = std::min(ceiling, sum_left - 5 * vertices_left);
    for (int value = upper; value >= 5; --value) {
        if (sum_left - value > value * vertices_left) continue;
        d[pos] = value;
        degree_sequences_rec(pos + 1, value, sum_left - value, d, sequences);
    }
}

std::vector<Degrees> degree_sequences() {
    std::vector<Degrees> result;
    Degrees d{};
    degree_sequences_rec(0, 2 * M - 30, 2 * M, d, result);
    return result;
}

using Graph = std::array<uint8_t, 21>;
Graph canonical(const Graph &g);
Degrees graph_degrees(const Graph &g) {
    Degrees d{};
    for (int i = 0; i < 21; ++i) {
        d[pairs[i].u] += g[i];
        d[pairs[i].v] += g[i];
    }
    return d;
}
Graph relabel_graph(const Graph &g, const std::array<int, N> &new_to_old) {
    Graph h{};
    for (int i = 0; i < 21; ++i) {
        int u = new_to_old[pairs[i].u], v = new_to_old[pairs[i].v];
        if (u > v) std::swap(u, v);
        int j = 0;
        while (pairs[j].u != u || pairs[j].v != v) ++j;
        h[i] = g[j];
    }
    return h;
}
Graph canonical_unfixed(const Graph &g) {
    Degrees d = graph_degrees(g);
    std::array<int, N> order{0, 1, 2, 3, 4, 5, 6};
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return d[a] > d[b]; });
    Graph sorted = relabel_graph(g, order);
    build_permutations(graph_degrees(sorted));
    return canonical(sorted);
}
struct Enumeration {
    uint64_t realizations = 0, labelled = 0;
    std::vector<Graph> classes;
};

bool edge_connected_5(const Graph &g) {
    for (int x : cuts) {
        int crossing = 0;
        for (int i = 0; i < 21; ++i)
            if (((x >> pairs[i].u) ^ (x >> pairs[i].v)) & 1) crossing += g[i];
        if (crossing < 5) return false;
    }
    return true;
}

Graph canonical(const Graph &g) {
    Graph best = g;
    for (const auto &pi : perm_indices) {
        int i = 0;
        for (; i < 21; ++i) {
            uint8_t a = g[pi[i]], b = best[i];
            if (a < b) break;
            if (a > b) { i = 21; break; }
        }
        if (i < 21) for (int j = 0; j < 21; ++j) best[j] = g[pi[j]];
    }
    return best;
}

void enumerate_rec(int pos, Graph &g, std::array<int, N> &rem, Enumeration &out,
                   std::map<Graph, bool> &unique) {
    if (pos == 21) {
        if (std::any_of(rem.begin(), rem.end(), [](int x) { return x != 0; })) return;
        ++out.realizations;
        if (!edge_connected_5(g)) return;
        ++out.labelled;
        unique.emplace(canonical(g), true);
        return;
    }
    int u = pairs[pos].u, v = pairs[pos].v;
    int max_m = std::min(rem[u], rem[v]);
    for (int m = 0; m <= max_m; ++m) {
        g[pos] = static_cast<uint8_t>(m);
        rem[u] -= m; rem[v] -= m;
        bool possible = true;
        for (int z = 0; z < N && possible; ++z) {
            int capacity = 0;
            for (int j = pos + 1; j < 21; ++j)
                if (pairs[j].u == z || pairs[j].v == z)
                    capacity += rem[pairs[j].u == z ? pairs[j].v : pairs[j].u];
            if (rem[z] > capacity) possible = false;
        }
        if (possible) enumerate_rec(pos + 1, g, rem, out, unique);
        rem[u] += m; rem[v] += m;
    }
}

Enumeration enumerate_graphs() {
    Enumeration out;
    Graph g{};
    std::map<Graph, bool> unique;
    auto sequences = degree_sequences();
    for (const Degrees &d : sequences) {
        build_permutations(d);
        Degrees rem = d;
        enumerate_rec(0, g, rem, out, unique);
    }
    for (const auto &[key, _] : unique) out.classes.push_back(key);
    std::cout << "M=" << M << " degree sequences: " << sequences.size() << '\n'
              << "degree-assigned realizations before lambda filter: " << out.realizations << '\n'
              << "degree-assigned lambda>=5 graphs: " << out.labelled << '\n'
              << "isomorphism classes: " << out.classes.size() << std::endl;
    if (M == 18 && (out.realizations != 747900 || out.labelled != 214430 || out.classes.size() != 372))
        throw std::runtime_error("graph enumeration checksum mismatch; orientation search stopped");
    return out;
}

struct CutInfo { U crossing = 0, zero_out = 0, one_out = 0; int total = 0; };
struct Engine {
    Graph graph;
    std::vector<Pair> edges;
    std::array<CutInfo, 63> ci{};
    std::vector<uint8_t> legal;
    std::vector<U> stamp, queue;
    U generation = 0;
    uint64_t legal_count = 0;

    explicit Engine(const Graph &g) : graph(g), edges(M), legal(S, 0), stamp(S, 0) {
        int e = 0;
        for (int i = 0; i < 21; ++i) for (int t = 0; t < g[i]; ++t) {
            if (e >= M) throw std::runtime_error("too many physical edges");
            edges[e++] = pairs[i];
        }
        if (e != M) throw std::runtime_error("physical edge count mismatch");
        for (int k = 0; k < 63; ++k) {
            int x = cuts[k];
            for (int j = 0; j < M; ++j) {
                int u = edges[j].u, v = edges[j].v;
                if (((x >> u) ^ (x >> v)) & 1) {
                    ci[k].crossing |= U(1) << j;
                    if (x & (1 << u)) ci[k].zero_out |= U(1) << j;
                    else ci[k].one_out |= U(1) << j;
                    ++ci[k].total;
                }
            }
        }
    }

    bool check(U state) const {
        for (const auto &c : ci) {
            int out = std::popcount((~state) & c.zero_out)
                    + std::popcount(state & c.one_out);
            if (out < 2 || c.total - out < 2) return false;
        }
        return true;
    }
    bool slow_check(U state) const {
        for (int x = 1; x < 127; ++x) {
            int out = 0;
            for (int e = 0; e < M; ++e) {
                int u = edges[e].u, v = edges[e].v;
                if (state & (U(1) << e)) std::swap(u, v);
                out += ((x >> u) & 1) && !((x >> v) & 1);
            }
            if (out < 2) return false;
        }
        return true;
    }
    void prepare(unsigned workers = 1) {
        legal_count = 0;
        workers = std::max(1u, std::min<unsigned>(workers, S));
        if (workers == 1) {
            for (U s = 0; s < S; ++s) {
                legal[s] = check(s);
                legal_count += legal[s];
            }
        } else {
            std::vector<uint64_t> counts(workers, 0);
            std::vector<std::thread> threads;
            threads.reserve(workers);
            U chunk = (S + workers - 1) / workers;
            try {
                for (unsigned id = 0; id < workers; ++id) threads.emplace_back([&, id] {
                    U begin = id * chunk, end = std::min<U>(begin + chunk, S);
                    uint64_t count = 0;
                    for (U s = begin; s < end; ++s) {
                        legal[s] = check(s);
                        count += legal[s];
                    }
                    counts[id] = count;
                });
            } catch (...) {
                for (auto &t : threads) t.join();
                throw;
            }
            for (auto &t : threads) t.join();
            for (uint64_t count : counts) legal_count += count;
        }
        queue.reserve(static_cast<size_t>(legal_count));
    }

    // Returns -1 only if the goal is unreachable in the legal flip graph.
    int distance(U source, U target, bool path, std::vector<U> *states = nullptr) {
        if (++generation == 0) { std::fill(stamp.begin(), stamp.end(), 0); generation = 1; }
        queue.clear(); queue.push_back(source); stamp[source] = generation;
        std::vector<U> pred;
        if (path) { pred.assign(S, 0); pred[source] = source; }
        size_t head = 0, layer_end = 1;
        int depth = 0;
        while (head < queue.size()) {
            while (head < layer_end) {
                U x = queue[head++];
                if (std::popcount(x ^ target) <= 2) {
                    if (path && states) {
                        states->clear();
                        for (U y = x;; y = pred[y]) {
                            states->push_back(y);
                            if (y == source) break;
                        }
                        std::reverse(states->begin(), states->end());
                    }
                    return depth;
                }
                for (int e = 0; e < M; ++e) {
                    U y = x ^ (U(1) << e);
                    if (legal[y] && stamp[y] != generation) {
                        stamp[y] = generation;
                        if (path) pred[y] = x;
                        queue.push_back(y);
                    }
                }
            }
            layer_end = queue.size();
            ++depth;
        }
        return -1;
    }
};

std::string hist_text(const Counts &h) {
    std::ostringstream os;
    bool first = true;
    for (auto [d, n] : h) { if (!first) os << ','; first = false; os << d << ':' << n; }
    return os.str();
}
std::string report_hist(const Counts &h) {
    std::ostringstream os;
    int max_d = h.empty() ? 6 : std::max(6, h.rbegin()->first);
    for (int d = 1; d <= max_d; ++d) {
        if (d > 1) os << ',';
        auto it = h.find(d);
        os << d << ':' << (it == h.end() ? 0 : it->second);
    }
    return os.str();
}
Counts parse_hist(const std::string &s) {
    Counts h;
    std::istringstream in(s);
    std::string part;
    while (std::getline(in, part, ',')) {
        if (part.empty()) continue;
        auto col = part.find(':');
        if (col == std::string::npos) throw std::runtime_error("bad histogram checkpoint");
        h[std::stoi(part.substr(0, col))] = std::stoull(part.substr(col + 1));
    }
    return h;
}

struct Progress {
    int next = 0, record = (M == 18 ? 5 : 0), edge_count = M;
    uint64_t graph_classes = 0;
    uint64_t legal_total = 0, pair_total = 0, records_bytes = 0;
    Counts hist, saved;
};
std::string serialize(const Progress &p) {
    std::ostringstream s;
    s << "version=2\nM=" << p.edge_count << "\ngraph_classes=" << p.graph_classes
      << "\nnext=" << p.next << "\nrecord=" << p.record
      << "\nlegal_total=" << p.legal_total << "\npair_total=" << p.pair_total
      << "\nrecords_bytes=" << p.records_bytes << "\nhist=" << hist_text(p.hist)
      << "\nsaved=" << hist_text(p.saved) << '\n';
    return s.str();
}
Progress deserialize(const std::string &text) {
    Progress p;
    std::map<std::string, std::string> kv;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        auto eq = line.find('=');
        if (eq != std::string::npos) kv[line.substr(0, eq)] = line.substr(eq + 1);
    }
    if (kv.at("version") != "1" && kv.at("version") != "2")
        throw std::runtime_error("unsupported checkpoint version");
    if (kv.at("version") == "1" && M != 18)
        throw std::runtime_error("legacy checkpoint belongs to M=18");
    if (kv.at("version") == "2") {
        p.edge_count = std::stoi(kv.at("M"));
        p.graph_classes = std::stoull(kv.at("graph_classes"));
        if (p.edge_count != M) throw std::runtime_error("checkpoint M mismatch");
    }
    p.next = std::stoi(kv.at("next")); p.record = std::stoi(kv.at("record"));
    p.legal_total = std::stoull(kv.at("legal_total"));
    p.pair_total = std::stoull(kv.at("pair_total"));
    p.records_bytes = std::stoull(kv.at("records_bytes"));
    p.hist = parse_hist(kv.at("hist")); p.saved = parse_hist(kv.at("saved"));
    return p;
}
void write_atomic(const fs::path &path, const std::string &contents) {
    fs::path temp = path; temp += ".tmp";
    { std::ofstream o(temp, std::ios::binary | std::ios::trunc); if (!o) throw std::runtime_error("cannot write " + temp.string()); o << contents; o.flush(); if (!o) throw std::runtime_error("write failed " + temp.string()); }
#ifdef _WIN32
    DWORD error = ERROR_SUCCESS;
    for (int attempt = 0; attempt < 20; ++attempt) {
        if (MoveFileExW(temp.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
            return;
        error = GetLastError();
        if (error != ERROR_ACCESS_DENIED && error != ERROR_SHARING_VIOLATION
            && error != ERROR_LOCK_VIOLATION) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(25 * (attempt + 1)));
    }
    throw std::runtime_error("atomic replace failed " + path.string()
                             + " Win32_error=" + std::to_string(error));
#else
    fs::rename(temp, path);
#endif
}

#ifdef _WIN32
struct SearchLock {
    HANDLE handle = INVALID_HANDLE_VALUE;
    explicit SearchLock(const fs::path &path) {
        handle = CreateFileW(path.c_str(), GENERIC_READ | GENERIC_WRITE, 0, nullptr,
                             OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle == INVALID_HANDLE_VALUE)
            throw std::runtime_error("cannot acquire search lock " + path.string()
                                     + " Win32_error=" + std::to_string(GetLastError())
                                     + "; another run may already be active");
    }
    SearchLock(const SearchLock &) = delete;
    SearchLock &operator=(const SearchLock &) = delete;
    ~SearchLock() { if (handle != INVALID_HANDLE_VALUE) CloseHandle(handle); }
};
#else
struct SearchLock { explicit SearchLock(const fs::path &) {} };
#endif
std::string read_all(const fs::path &p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("cannot read " + p.string());
    return std::string(std::istreambuf_iterator<char>(f), {});
}

void repair_files(const Progress &p) {
    fs::path rec = OUTPUT_DIR / "records.jsonl", csv = OUTPUT_DIR / "per_graph.csv";
    if (!fs::exists(rec)) { if (p.records_bytes) throw std::runtime_error("records file missing"); std::ofstream(rec); }
    if (fs::file_size(rec) < p.records_bytes) throw std::runtime_error("records shorter than checkpoint");
    fs::resize_file(rec, p.records_bytes);
    std::string rows = read_all(csv);
    size_t end = 0;
    for (int i = 0; i <= p.next; ++i) {
        end = rows.find('\n', end);
        if (end == std::string::npos) throw std::runtime_error("per_graph.csv shorter than checkpoint");
        ++end;
    }
    fs::resize_file(csv, end);
    std::string last;
    std::ifstream in(rec);
    std::string line;
    std::string needle = "\"distance\":" + std::to_string(p.record) + ",";
    while (std::getline(in, line)) if (line.find(needle) != std::string::npos) last = line;
    write_atomic(OUTPUT_DIR / "best_witness.txt", last.empty() ? "No committed record witness yet.\n" : last + "\n");
}

std::string witness_json(int gi, const Engine &eng, U d, U target, const std::array<int, 3> &mis,
                         int distance, const std::vector<U> &path) {
    std::ostringstream s;
    s << "{\"M\":" << M << ",\"graph_class_index\":" << gi << ",\"degrees\":[";
    for (int v = 0; v < N; ++v) {
        if (v) s << ',';
        int degree = 0;
        for (int i = 0; i < 21; ++i)
            if (pairs[i].u == v || pairs[i].v == v) degree += eng.graph[i];
        s << degree;
    }
    s << "],\"multiplicities\":[";
    for (int i = 0; i < 21; ++i) { if (i) s << ','; s << int(eng.graph[i]); }
    s << "],\"multiplicity_pair_endpoints\":[";
    for (int i = 0; i < 21; ++i) {
        if (i) s << ',';
        s << '[' << pairs[i].u << ',' << pairs[i].v << ']';
    }
    s << "],\"physical_edges\":[";
    for (int e = 0; e < M; ++e) { if (e) s << ','; s << '[' << eng.edges[e].u << ',' << eng.edges[e].v << ']'; }
    s << "],\"source\":" << d << ",\"target\":" << target << ",\"initial_mismatch_edges\":["
      << mis[0] << ',' << mis[1] << ',' << mis[2] << "],\"source_legal_neighbors\":[";
    bool first = true;
    for (int e = 0; e < M; ++e) if (eng.legal[d ^ (U(1) << e)]) {
        if (!first) s << ','; first = false;
        s << "{\"edge\":" << e << ",\"state\":" << (d ^ (U(1) << e)) << '}';
    }
    s << "],\"distance\":" << distance << ",\"path\":[";
    for (size_t i = 0; i < path.size(); ++i) {
        if (i) s << ',';
        int edge = -1;
        if (i) edge = std::countr_zero(path[i] ^ path[i - 1]);
        s << "{\"state\":" << path[i] << ",\"flipped_edge\":" << edge
          << ",\"mismatch\":" << std::popcount(path[i] ^ target)
          << ",\"move\":\"" << (i == 0 ? "START" : ((path[i - 1] ^ target) & (U(1) << edge) ? "CORRECTION" : "COMMON"))
          << "\",\"legal\":" << (eng.check(path[i]) ? "true" : "false") << '}';
    }
    s << "],\"final_mismatch_at_most_2\":" << (std::popcount(path.back() ^ target) <= 2 ? "true" : "false") << '}';
    return s.str();
}

void summary(const Progress &p, bool complete) {
    std::ostringstream s;
    int exact_max = p.hist.empty() ? 0 : p.hist.rbegin()->first;
    s << "COMPLETE_EXHAUSTIVE_M" << M << " = " << (complete ? "YES" : "NO") << '\n'
      << "graph_classes = " << (p.graph_classes ? std::to_string(p.graph_classes) : "PENDING_ENUMERATION") << '\n'
      << "graph_classes_processed = " << p.next << '/'
      << (p.graph_classes ? std::to_string(p.graph_classes) : "PENDING_ENUMERATION") << '\n'
      << "total_legal_states = " << p.legal_total << '\n'
      << "total_qualifying_ordered_pairs = " << p.pair_total << '\n'
      << "global_distance_histogram = " << report_hist(p.hist) << '\n'
      << "global_maximum_progress_distance = " << exact_max << '\n'
      << "witnesses_attaining_maximum = " << (p.hist.count(exact_max) ? p.hist.at(exact_max) : 0) << '\n';
    write_atomic(OUTPUT_DIR / "summary.txt", s.str());
}

void self_test(const Graph &g, uint64_t class_count) {
    Engine eng(g); eng.prepare();
    if (!eng.legal_count) throw std::runtime_error("self-test found no legal states");
    if (M == 18) {
        for (U s = 0; s < S; ++s) if (eng.legal[s] != eng.slow_check(s))
            throw std::runtime_error("legal-state self-test mismatch");
    } else {
        for (U s = 0; s < std::min<U>(S, 4096); ++s)
            if (eng.legal[s] != eng.slow_check(s)) throw std::runtime_error("legal-state self-test mismatch");
        for (U s = S - 4096; s < S; ++s)
            if (eng.legal[s] != eng.slow_check(s)) throw std::runtime_error("legal-state self-test mismatch");
        int checked_legal = 0;
        for (U s = 0; s < S && checked_legal < 256; ++s) if (eng.legal[s]) {
            if (!eng.slow_check(s)) throw std::runtime_error("legal-state self-test mismatch");
            ++checked_legal;
        }
    }
    U d = 0, target = 0;
    bool found = false;
    for (U s = 0; s < S && !found; ++s) if (eng.legal[s]) {
        for (int a = 0; a < M && !found; ++a) for (int b = a + 1; b < M && !found; ++b)
            for (int c = b + 1; c < M && !found; ++c) {
                U t = s ^ (U(1) << a) ^ (U(1) << b) ^ (U(1) << c);
                if (eng.legal[t] && (M != 18 || (!eng.legal[s ^ (U(1) << a)]
                    && !eng.legal[s ^ (U(1) << b)] && !eng.legal[s ^ (U(1) << c)])))
                    { d = s; target = t; found = true; }
            }
    }
    if (!found) throw std::runtime_error("self-test graph lacks legal Hamming-3 pair");
    std::vector<U> path;
    int dist = eng.distance(d, target, true, &path);
    if (dist < 0 || int(path.size()) != dist + 1 || path.front() != d
        || std::popcount(path.back() ^ target) > 2) throw std::runtime_error("BFS path self-test failure");
    for (size_t i = 1; i < path.size(); ++i)
        if (std::popcount(path[i] ^ path[i - 1]) != 1 || !eng.legal[path[i]]
            || (i + 1 < path.size() && std::popcount(path[i] ^ target) <= 2))
            throw std::runtime_error("BFS shortest-path self-test failure");
    // Independent distance-only reference uses ordinary visited bytes and layer BFS.
    std::vector<uint8_t> seen(S, 0);
    std::vector<U> q{d}; seen[d] = 1;
    size_t head = 0, end = 1; int ref = 0;
    bool hit = false;
    while (head < q.size() && !hit) {
        while (head < end && !hit) {
            U x = q[head++];
            if (std::popcount(x ^ target) <= 2) { hit = true; break; }
            for (int e = 0; e < M; ++e) {
                U y = x ^ (U(1) << e);
                if (eng.legal[y] && !seen[y]) { seen[y] = 1; q.push_back(y); }
            }
        }
        if (!hit) { end = q.size(); ++ref; }
    }
    if (!hit || ref != dist) throw std::runtime_error("BFS reference distance mismatch");
    Progress p; p.next = 17; p.graph_classes = class_count;
    p.legal_total = 1234; p.pair_total = 567;
    p.hist[3] = 10; p.hist[5] = 2; p.saved[5] = 1; p.records_bytes = 99;
    fs::create_directories(OUTPUT_DIR);
    fs::path test = OUTPUT_DIR / "self_test_checkpoint.tmp";
    write_atomic(test, serialize(p));
    Progress r = deserialize(read_all(test)); fs::remove(test);
    if (serialize(p) != serialize(r)) throw std::runtime_error("checkpoint serialization mismatch");
#ifdef _WIN32
    fs::path retry_test = OUTPUT_DIR / "self_test_atomic_replace.tmp";
    write_atomic(retry_test, "old\n");
    HANDLE held = CreateFileW(retry_test.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                              OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (held == INVALID_HANDLE_VALUE) throw std::runtime_error("atomic replace self-test lock failed");
    std::thread release([held] {
        std::this_thread::sleep_for(std::chrono::milliseconds(150));
        CloseHandle(held);
    });
    try { write_atomic(retry_test, "new\n"); }
    catch (...) { release.join(); throw; }
    release.join();
    if (read_all(retry_test) != "new\n") throw std::runtime_error("atomic replace retry self-test failed");
    fs::remove(retry_test);
#endif
    if (!fs::exists(OUTPUT_DIR / "checkpoint.txt")) {
        Progress initial; initial.graph_classes = class_count;
        write_atomic(OUTPUT_DIR / "records.jsonl", "");
        write_atomic(OUTPUT_DIR / "per_graph.csv", "graph_class_index,legal_orientation_count,qualifying_ordered_pair_count,maximum_distance,distance_histogram\n");
        write_atomic(OUTPUT_DIR / "best_witness.txt", "No record witness yet.\n");
        write_atomic(OUTPUT_DIR / "checkpoint.txt", serialize(initial));
        summary(initial, false);
    } else if (class_count == 0) {
        Progress existing = deserialize(read_all(OUTPUT_DIR / "checkpoint.txt"));
        if (existing.next == 0 && existing.graph_classes == 0) summary(existing, false);
    }
    std::cout << "SELF_TEST_OK M=" << M << " legal_states=" << eng.legal_count
              << " BFS_distance=" << dist << std::endl;
}

// One source orientation is an independent unit of exact pair search.
template <class OnPair>
void scan_source(Engine &eng, U d, OnPair &&on_pair) {
    if (!eng.legal[d]) return;
    U blocked = 0;
    for (int e = 0; e < M; ++e)
        if (!eng.legal[d ^ (U(1) << e)]) blocked |= U(1) << e;
    if (std::popcount(blocked) < 3) return;
    for (int a = 0; a < M - 2; ++a) if (blocked & (U(1) << a))
        for (int b = a + 1; b < M - 1; ++b) if (blocked & (U(1) << b))
            for (int c = b + 1; c < M; ++c) if (blocked & (U(1) << c)) {
                U target = d ^ (U(1) << a) ^ (U(1) << b) ^ (U(1) << c);
                if (!eng.legal[target]) continue;
                int dist = eng.distance(d, target, false);
                if (dist < 0) throw std::runtime_error("qualifying pair has unreachable progress goal");
                on_pair(d, target, std::array<int, 3>{a, b, c}, dist);
            }
}

struct ClassStats { Counts hist; uint64_t pairs = 0; };

unsigned m19_worker_count() {
    unsigned logical = std::thread::hardware_concurrency();
    if (!logical) logical = 1;
    unsigned reserve = logical >= 6 ? 2 : (logical >= 2 ? 1 : 0);
    return std::max(1u, logical - reserve);
}


ClassStats scan_class_parallel19(const Engine &eng, int gi, Progress &p,
                                 std::ofstream &record_out, unsigned requested_workers,
                                 const fs::path &best_witness_path,
                                 const std::vector<U> *source_override = nullptr) {
    std::vector<U> sources;
    if (source_override) sources = *source_override;
    else {
        sources.reserve(static_cast<size_t>(eng.legal_count));
        for (U d = 0; d < S; ++d) if (eng.legal[d]) sources.push_back(d);
    }
    if (sources.empty()) return {};
    unsigned workers = std::min<unsigned>(requested_workers, static_cast<unsigned>(sources.size()));
    std::vector<ClassStats> parts(workers);
    std::atomic<size_t> next_source{0};
    std::atomic<bool> abort{false};
    std::atomic<int> live_record{p.record};
    Counts live_saved = p.saved;
    std::mutex record_mutex, failure_mutex;
    std::exception_ptr failure;

    auto run_worker = [&](unsigned id) {
        try {
            Engine worker = eng; // Private generation stamps and BFS queue per worker.
            worker.queue.reserve(static_cast<size_t>(worker.legal_count));
            while (!abort.load(std::memory_order_relaxed)) {
                size_t begin = next_source.fetch_add(8, std::memory_order_relaxed);
                if (begin >= sources.size()) break;
                size_t end = std::min(begin + 8, sources.size());
                for (size_t i = begin; i < end && !abort.load(std::memory_order_relaxed); ++i) {
                    scan_source(worker, sources[i], [&](U d, U target,
                                                        const std::array<int, 3> &mis, int dist) {
                        ++parts[id].pairs;
                        ++parts[id].hist[dist];
                        if (dist < live_record.load(std::memory_order_relaxed)) return;
                        std::lock_guard<std::mutex> lock(record_mutex);
                        int current = live_record.load(std::memory_order_relaxed);
                        if (dist < current) return;
                        bool new_record = dist > current;
                        if (new_record) {
                            live_record.store(dist, std::memory_order_relaxed);
                            if (gi >= 0) std::cout << "NEW_RECORD graph=" << gi << " source=" << d
                                                   << " target=" << target << " distance=" << dist << std::endl;
                        }
                        if (live_saved[dist] >= 100) return;
                        std::vector<U> path;
                        int again = worker.distance(d, target, true, &path);
                        if (again != dist) throw std::runtime_error("witness BFS mismatch");
                        std::string w = witness_json(gi, worker, d, target, mis, dist, path);
                        record_out << w << '\n'; record_out.flush();
                        if (!record_out) throw std::runtime_error("record flush failed");
                        ++live_saved[dist];
                        if (new_record) {
                            write_atomic(best_witness_path, w + "\n");
                            if (gi >= 0) std::cout << w << std::endl;
                        }
                    });
                }
            }
        } catch (...) {
            abort.store(true, std::memory_order_relaxed);
            std::lock_guard<std::mutex> lock(failure_mutex);
            if (!failure) failure = std::current_exception();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(workers);
    try {
        for (unsigned id = 0; id < workers; ++id) threads.emplace_back(run_worker, id);
    } catch (...) {
        abort.store(true, std::memory_order_relaxed);
        for (auto &t : threads) t.join();
        throw;
    }
    for (auto &t : threads) t.join();
    if (failure) std::rethrow_exception(failure);
    ClassStats total;
    for (const auto &part : parts) {
        total.pairs += part.pairs;
        for (auto [dist, count] : part.hist) total.hist[dist] += count;
    }
    p.record = live_record.load(std::memory_order_relaxed);
    p.saved = std::move(live_saved);
    return total;
}

void search(const Enumeration &en, bool fresh) {
    fs::create_directories(OUTPUT_DIR);
    SearchLock search_lock(OUTPUT_DIR / "search.lock");
    Progress p; p.graph_classes = en.classes.size();
    fs::path cp = OUTPUT_DIR / "checkpoint.txt", rec = OUTPUT_DIR / "records.jsonl", csv = OUTPUT_DIR / "per_graph.csv";
    if (fresh || !fs::exists(cp)) {
        write_atomic(rec, "");
        write_atomic(csv, "graph_class_index,legal_orientation_count,qualifying_ordered_pair_count,maximum_distance,distance_histogram\n");
        write_atomic(OUTPUT_DIR / "best_witness.txt", "No record witness yet.\n");
        write_atomic(cp, serialize(p));
    } else {
        p = deserialize(read_all(cp));
        if (p.graph_classes && p.graph_classes != en.classes.size())
            throw std::runtime_error("checkpoint graph class count mismatch");
        p.graph_classes = en.classes.size();
        if (p.next < 0 || static_cast<size_t>(p.next) > en.classes.size())
            throw std::runtime_error("bad checkpoint graph index");
        repair_files(p);
    }
    summary(p, static_cast<size_t>(p.next) == en.classes.size());
    std::ofstream record_out(rec, std::ios::binary | std::ios::app);
    if (!record_out) throw std::runtime_error("cannot open records.jsonl");
    unsigned workers = M == 19 ? m19_worker_count() : 1;
    if (M >= 20) workers = m19_worker_count();
    if (M == 19) std::cout << "M=19 orientation search workers=" << workers
                           << " (logical processors=" << std::thread::hardware_concurrency() << ")\n";
    else if (M >= 20) std::cout << "M=" << M << " orientation search workers=" << workers
                                << " (logical processors=" << std::thread::hardware_concurrency() << ")\n";
    for (int gi = p.next; static_cast<size_t>(gi) < en.classes.size(); ++gi) {
        Engine eng(en.classes[gi]); eng.prepare(workers);
        ClassStats class_stats;
        if (M == 19 && workers > 1) {
            class_stats = scan_class_parallel19(eng, gi, p, record_out, workers,
                                                OUTPUT_DIR / "best_witness.txt");
        } else if (M >= 20 && workers > 1) {
            class_stats = scan_class_parallel19(eng, gi, p, record_out, workers,
                                                OUTPUT_DIR / "best_witness.txt");
        } else {
            for (U d = 0; d < S; ++d) {
                scan_source(eng, d, [&](U source, U target,
                                       const std::array<int, 3> &mis, int dist) {
                    ++class_stats.pairs; ++class_stats.hist[dist];
                    if (dist > p.record) {
                        p.record = dist;
                        std::cout << "NEW_RECORD graph=" << gi << " source=" << source
                                  << " target=" << target << " distance=" << dist << std::endl;
                    }
                    if (dist >= p.record && p.saved[dist] < 100) {
                        std::vector<U> path;
                        int again = eng.distance(source, target, true, &path);
                        if (again != dist) throw std::runtime_error("witness BFS mismatch");
                        std::string w = witness_json(gi, eng, source, target, mis, dist, path);
                        record_out << w << '\n'; record_out.flush();
                        if (!record_out) throw std::runtime_error("record flush failed");
                        ++p.saved[dist];
                        if (p.saved[dist] == 1) {
                            write_atomic(OUTPUT_DIR / "best_witness.txt", w + "\n");
                            if (dist > 5) std::cout << w << std::endl;
                        }
                    }
                });
            }
        }
        for (auto [dist, count] : class_stats.hist) p.hist[dist] += count;
        p.pair_total += class_stats.pairs;
        p.legal_total += eng.legal_count;
        int local_max = class_stats.hist.empty() ? 0 : class_stats.hist.rbegin()->first;
        { std::ofstream f(csv, std::ios::app); if (!f) throw std::runtime_error("cannot append per_graph.csv");
          f << gi << ',' << eng.legal_count << ',' << class_stats.pairs << ',' << local_max
            << ",\"" << report_hist(class_stats.hist) << "\"\n"; f.flush(); if (!f) throw std::runtime_error("per_graph flush failed"); }
        record_out.flush();
        p.records_bytes = fs::file_size(rec);
        p.next = gi + 1;
        write_atomic(cp, serialize(p));
        summary(p, static_cast<size_t>(p.next) == en.classes.size());
        std::cout << "completed graph class " << gi << '/' << (en.classes.size() - 1) << "; legal=" << eng.legal_count
                  << " pairs=" << class_stats.pairs << " max=" << local_max << std::endl;
    }
}

void self_test_parallel19(const Graph &g) {
    Engine eng(g); eng.prepare(m19_worker_count());
    Engine serial_precompute(g); serial_precompute.prepare();
    if (eng.legal_count != serial_precompute.legal_count || eng.legal != serial_precompute.legal)
        throw std::runtime_error("parallel legal-state precompute differs from serial self-test");
    std::vector<U> sources;
    int legal_examined = 0;
    for (U d = 0; d < S && sources.size() < 16 && legal_examined < 20000; ++d) {
        if (!eng.legal[d]) continue;
        ++legal_examined;
        U blocked = 0;
        for (int e = 0; e < M; ++e)
            if (!eng.legal[d ^ (U(1) << e)]) blocked |= U(1) << e;
        if (std::popcount(blocked) < 3) continue;
        bool qualifies = false;
        for (int a = 0; a < M - 2 && !qualifies; ++a) if (blocked & (U(1) << a))
            for (int b = a + 1; b < M - 1 && !qualifies; ++b) if (blocked & (U(1) << b))
                for (int c = b + 1; c < M; ++c) if (blocked & (U(1) << c))
                    if (eng.legal[d ^ (U(1) << a) ^ (U(1) << b) ^ (U(1) << c)]) {
                        qualifies = true; break;
                    }
        if (qualifies) sources.push_back(d);
    }
    if (sources.empty()) throw std::runtime_error("parallel self-test found no qualifying M=19 sources");
    Engine serial = eng;
    serial.queue.reserve(static_cast<size_t>(serial.legal_count));
    ClassStats expected;
    for (U d : sources) scan_source(serial, d, [&](U, U, const std::array<int, 3> &, int dist) {
        ++expected.pairs; ++expected.hist[dist];
    });
    fs::path temp_records = OUTPUT_DIR / "parallel_self_test_records.tmp";
    fs::path temp_best = OUTPUT_DIR / "parallel_self_test_best.tmp";
    Progress test_progress;
    ClassStats actual;
    {
        std::ofstream f(temp_records, std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("parallel self-test record file open failed");
        actual = scan_class_parallel19(eng, -1, test_progress, f,
                                       m19_worker_count(), temp_best, &sources);
    }
    bool records_exist = fs::exists(temp_records) && fs::file_size(temp_records) > 0;
    bool best_exists = fs::exists(temp_best) && fs::file_size(temp_best) > 0;
    fs::remove(temp_records);
    if (fs::exists(temp_best)) fs::remove(temp_best);
    if (expected.pairs == 0 || expected.pairs != actual.pairs
        || expected.hist != actual.hist || !records_exist || !best_exists
        || test_progress.record != expected.hist.rbegin()->first)
        throw std::runtime_error("parallel M=19 search differs from serial self-test");
    std::cout << "PARALLEL_SELF_TEST_OK M=19 sources=" << sources.size()
              << " qualifying_pairs=" << actual.pairs
              << " workers=" << std::min<unsigned>(m19_worker_count(), sources.size()) << std::endl;
}

void self_test_parallel20_22(const Graph &g) {
    Engine eng(g); eng.prepare(m19_worker_count());
    Engine serial_precompute(g); serial_precompute.prepare();
    if (eng.legal_count != serial_precompute.legal_count || eng.legal != serial_precompute.legal)
        throw std::runtime_error("parallel legal-state precompute differs from serial self-test");
    std::vector<U> sources;
    int legal_examined = 0;
    for (U d = 0; d < S && sources.size() < 32 && legal_examined < 20000; ++d) {
        if (!eng.legal[d]) continue;
        ++legal_examined;
        U blocked = 0;
        for (int e = 0; e < M; ++e)
            if (!eng.legal[d ^ (U(1) << e)]) blocked |= U(1) << e;
        if (std::popcount(blocked) < 3) continue;
        bool qualifies = false;
        for (int a = 0; a < M - 2 && !qualifies; ++a) if (blocked & (U(1) << a))
            for (int b = a + 1; b < M - 1 && !qualifies; ++b) if (blocked & (U(1) << b))
                for (int c = b + 1; c < M; ++c) if (blocked & (U(1) << c))
                    if (eng.legal[d ^ (U(1) << a) ^ (U(1) << b) ^ (U(1) << c)]) {
                        qualifies = true; break;
                    }
        if (qualifies) sources.push_back(d);
    }
    if (sources.empty()) throw std::runtime_error("parallel self-test found no qualifying M>=20 sources");
    Engine serial = eng;
    serial.queue.reserve(static_cast<size_t>(serial.legal_count));
    ClassStats expected;
    for (U d : sources) scan_source(serial, d, [&](U, U, const std::array<int, 3> &, int dist) {
        ++expected.pairs; ++expected.hist[dist];
    });
    fs::path temp_records = OUTPUT_DIR / "parallel_self_test_records.tmp";
    fs::path temp_best = OUTPUT_DIR / "parallel_self_test_best.tmp";
    Progress test_progress;
    ClassStats actual;
    {
        std::ofstream f(temp_records, std::ios::binary | std::ios::trunc);
        if (!f) throw std::runtime_error("parallel self-test record file open failed");
        actual = scan_class_parallel19(eng, -1, test_progress, f,
                                       m19_worker_count(), temp_best, &sources);
    }
    bool records_exist = fs::exists(temp_records) && fs::file_size(temp_records) > 0;
    bool best_exists = fs::exists(temp_best) && fs::file_size(temp_best) > 0;
    fs::remove(temp_records);
    if (fs::exists(temp_best)) fs::remove(temp_best);
    if (expected.pairs == 0 || expected.pairs != actual.pairs
        || expected.hist != actual.hist || !records_exist || !best_exists
        || test_progress.record != expected.hist.rbegin()->first)
        throw std::runtime_error("parallel M>=20 search differs from serial self-test");
    std::cout << "PARALLEL_SELF_TEST_OK M=" << M << " sources=" << sources.size()
              << " qualifying_pairs=" << actual.pairs
              << " workers=" << std::min<unsigned>(m19_worker_count(), sources.size()) << std::endl;
}

void self_test_selected_M() {
    if (M == 18) {
        Enumeration en = enumerate_graphs();
        self_test(en.classes.front(), en.classes.size());
        return;
    }
    int selected_M = M;
    M = 18; S = U(1) << M;
    Enumeration reference = enumerate_graphs();
    M = selected_M; S = U(1) << M;
    std::cout << "self-test M=" << M
              << " uses a synthetic graph; M=" << M << " class enumeration is reserved for the full run\n";
    Graph g = reference.classes.front();
    g[0] += static_cast<uint8_t>(M - 18);
    if (!edge_connected_5(g)) throw std::runtime_error("synthetic self-test graph failed lambda>=5");
    Degrees d = graph_degrees(g);
    if (std::accumulate(d.begin(), d.end(), 0) != 2 * M
        || !std::is_sorted(d.begin(), d.end(), std::greater<int>()))
        throw std::runtime_error("synthetic degree sequence invalid");
    auto sequences = degree_sequences();
    if (std::find(sequences.begin(), sequences.end(), d) == sequences.end())
        throw std::runtime_error("degree sequence enumeration missed synthetic graph");
    Graph base = canonical_unfixed(g);
    std::array<int, N> order{0, 1, 2, 3, 4, 5, 6};
    for (int shift = 0; shift < N; ++shift) {
        Graph relabelled = relabel_graph(g, order);
        if (canonical_unfixed(relabelled) != base)
            throw std::runtime_error("isomorphism canonicalization self-test failed");
        std::rotate(order.begin(), order.begin() + 1, order.end());
    }
    std::reverse(order.begin(), order.end());
    if (canonical_unfixed(relabel_graph(g, order)) != base)
        throw std::runtime_error("isomorphism reversal self-test failed");
    self_test(g, 0);
    if (M == 19) self_test_parallel19(g);
    else if (M >= 20) self_test_parallel20_22(g);
}

int main(int argc, char **argv) {
    try {
        bool fresh = false, test = false;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--fresh") fresh = true;
            else if (a == "--self-test") test = true;
            else if (a == "--M") {
                if (i + 1 >= argc) throw std::runtime_error("--M requires 18, 19, 20, 21, or 22");
                M = std::stoi(argv[++i]);
                if (M < 18 || M > 22) throw std::runtime_error("--M must be 18, 19, 20, 21, or 22");
            }
            else throw std::runtime_error("unknown argument: " + a);
        }
        S = U(1) << M;
        OUTPUT_DIR = fs::path(ROOT) / ("m" + std::to_string(M) + "_results");
        init_tables();
        if (test) self_test_selected_M();
        else {
            Enumeration en = enumerate_graphs();
            search(en, fresh);
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}

