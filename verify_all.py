#!/usr/bin/env python3
"""
Exact symbolic verification script accompanying:
  "A Spectral Dichotomy for Dilworth Multiplicity Bounds of the Signless Laplacian
   and Distance Signless Laplacian on C-Graphs and Cographs"

Verifies:
  1. Theorem 3.1: Infinite counterexample family G_k = C(k, 1, ..., 1) for k = 1..5
     directly on the full n x n integer matrix D^Q(G_k).
  2. Proposition 3.2: Exhaustive exact polynomial gcd(p, p') check of all 127
     connected C-graphs of order 2 <= n <= 8, proving H_8 = C(1^8) is the unique
     smallest counterexample to Mandal-Manna Conjecture 5.5.
  3. Theorem 5.4: Full n x n exact integer verification of B = K_36 V (3 * H_8)
     (n=60, nabla=3, mult_Q(40)=4) and F_s = K_1 V (s * B) for s = 2, 3
     (n=121, 181; mult_Q(41) >= 4s-1 and mult_{D^Q}(120s-43) >= 4s-1 = (4/3)nabla - 1).
  4. Table 2: Exhaustive exact symbolic verification of mult_Q <= nabla and
     mult_{D^Q} <= nabla + 1 across all 1,170 connected cographs of order 2 <= n <= 9.

Requirements: numpy, sympy (Python 3.8+)
"""

from itertools import combinations_with_replacement, product
import numpy as np
import sympy as sp


def build_cgraph_recursive(alphas):
    """Build adjacency matrix of C(alpha_1, ..., alpha_m) from the recursive definition."""
    n_cur = alphas[0]
    A = np.zeros((n_cur, n_cur), dtype=int)
    for a in alphas[1:]:
        Ka = np.ones((a, a), dtype=int) - np.eye(a, dtype=int)
        U = np.block([
            [A, np.zeros((n_cur, a), dtype=int)],
            [np.zeros((a, n_cur), dtype=int), Ka],
        ])
        n_new = n_cur + a
        A = np.ones((n_new, n_new), dtype=int) - np.eye(n_new, dtype=int) - U
        n_cur = n_new
    return A


def verify_theorem_3_1():
    print("[1/4] Verifying Theorem 3.1 (infinite family G_k, k = 1..5)...")
    for k in range(1, 6):
        m = 2 * k + 6
        alphas = [k] + [1] * (m - 1)
        A = build_cgraph_recursive(alphas)
        n = A.shape[0]
        assert n == 3 * k + 5
        D = 2 * (np.ones((n, n), dtype=int) - np.eye(n, dtype=int)) - A
        DQ = np.diag(np.sum(D, axis=1)) + D
        theta = 4 * k + 4
        x1 = np.array([-1] * k + [-1, 1, k] + [0] * (m - 4), dtype=int)
        x2 = np.array(
            [k + 1] * k + [k + 1, -(k + 1), 0] + [(-k if i % 2 == 0 else k) for i in range(m - 4)],
            dtype=int,
        )
        assert np.array_equal(DQ @ x1, theta * x1)
        assert np.array_equal(DQ @ x2, theta * x2)
        assert np.linalg.matrix_rank(np.column_stack([x1, x2]).astype(float)) == 2
    print("      PASSED: Exact integer eigenvector verification for k = 1..5.")


def compositions_even(n):
    res = []
    def bt(rem, cur):
        if rem == 0:
            if len(cur) >= 2 and len(cur) % 2 == 0:
                res.append(tuple(cur))
            return
        for a in range(1, rem + 1):
            bt(rem - a, cur + [a])
    bt(n, [])
    return res


def verify_proposition_3_2():
    print("[2/4] Verifying Proposition 3.2 (all 127 connected C-graphs, 2 <= n <= 8)...")
    x = sp.Symbol("x")
    total = 0
    violations = []
    for n in range(2, 9):
        for alphas in compositions_even(n):
            total += 1
            m = len(alphas)
            d = []
            for i in range(m):
                idx = i + 1
                deg = (alphas[i] - 1) if (idx % 2 == 1) else 0
                for j in range(m):
                    if j != i and max(idx, j + 1) % 2 == 0:
                        deg += alphas[j]
                d.append(deg)
            Tr = [2 * (n - 1) - deg for deg in d]
            TQ = {Tr[i] - (1 if (i + 1) % 2 == 1 else 2) for i in range(m) if alphas[i] >= 2}
            nabla = 1 + sum(1 for a in alphas[2:] if a >= 2)
            Q = sp.zeros(m, m)
            for i in range(m):
                idx = i + 1
                Q[i, i] = Tr[i] + (alphas[i] - 1) * (1 if idx % 2 == 1 else 2)
                for j in range(m):
                    if j != i:
                        Q[i, j] = alphas[j] * (1 if max(idx, j + 1) % 2 == 0 else 2)
            p = sp.Poly(Q.charpoly(x).as_expr(), x)
            if sp.gcd(p, p.diff()).degree() > 0:
                for e, mult in Q.eigenvals().items():
                    if e not in TQ and mult > nabla:
                        violations.append((alphas, nabla, mult, e))
    assert total == 127
    assert violations == [((1, 1, 1, 1, 1, 1, 1, 1), 1, 2, 8)]
    print(f"      PASSED: {total} C-graphs checked; unique violation is C(1^8) = H_8 (theta=8, mult=2 > nabla=1).")


def verify_theorem_5_4():
    print("[3/4] Verifying Theorem 5.4 (linear separation B = K_36 V 3*H_8 and F_s = K_1 V s*B)...")
    H8 = build_cgraph_recursive([1] * 8)
    w = 36
    nB = w + 24
    AB = np.zeros((nB, nB), dtype=int)
    AB[:w, :w] = 1 - np.eye(w, dtype=int)
    AB[:w, w:] = 1
    AB[w:, :w] = 1
    for t in range(3):
        AB[w + 8 * t : w + 8 * (t + 1), w + 8 * t : w + 8 * (t + 1)] = H8
    QB = np.diag(np.sum(AB, axis=1)) + AB
    x1 = np.array([-1, -1, 1, 1, 0, 0, 0, 0], dtype=int)
    x2 = np.array([2, 2, -2, 0, -1, 1, -1, 1], dtype=int)

    vecsB = []
    for t in range(3):
        v = np.zeros(nB, dtype=int)
        v[w + 8 * t : w + 8 * (t + 1)] = x1
        assert np.array_equal(QB @ v, 40 * v)
        vecsB.append(v)
    v4 = np.zeros(nB, dtype=int)
    v4[:w] = -1
    for t in range(3):
        v4[w + 8 * t : w + 8 * (t + 1)] = 9 * x2
    assert np.array_equal(QB @ v4, 40 * v4)
    vecsB.append(v4)
    assert np.linalg.matrix_rank(np.column_stack(vecsB).astype(float)) == 4

    for s in [2, 3]:
        nF = 1 + s * nB
        AF = np.zeros((nF, nF), dtype=int)
        AF[0, 1:] = 1
        AF[1:, 0] = 1
        for j in range(s):
            AF[1 + j * nB : 1 + (j + 1) * nB, 1 + j * nB : 1 + (j + 1) * nB] = AB
        QF = np.diag(np.sum(AF, axis=1)) + AF
        DF = 2 * (np.ones((nF, nF), dtype=int) - np.eye(nF, dtype=int)) - AF
        DQF = np.diag(np.sum(DF, axis=1)) + DF
        thetaF = 120 * s - 43
        vecsF = []
        for j in range(s):
            for t in range(3):
                v = np.zeros(nF, dtype=int)
                off = 1 + j * nB + w + 8 * t
                v[off : off + 8] = x1
                assert int(np.sum(v)) == 0
                assert np.array_equal(QF @ v, 41 * v)
                assert np.array_equal(DQF @ v, thetaF * v)
                vecsF.append(v)
        for j in range(1, s):
            v = np.zeros(nF, dtype=int)
            v[1 : 1 + nB] = v4
            v[1 + j * nB : 1 + (j + 1) * nB] = -v4
            assert int(np.sum(v)) == 0
            assert np.array_equal(QF @ v, 41 * v)
            assert np.array_equal(DQF @ v, thetaF * v)
            vecsF.append(v)
        assert np.linalg.matrix_rank(np.column_stack(vecsF).astype(float)) == 4 * s - 1
    print("      PASSED: B (n=60, nabla=3, mult_Q=4) and F_2, F_3 (mult_Q, mult_DQ >= (4/3)nabla - 1) verified.")


def partitions(n, max_part):
    if n == 0:
        yield ()
        return
    for p in range(min(n, max_part), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def gen_from_parts(pool, part_tuple):
    groups = []
    for p in part_tuple:
        if groups and groups[-1][0] == p:
            groups[-1][1] += 1
        else:
            groups.append([p, 1])
    lists = [list(combinations_with_replacement(pool[sz], cnt)) for sz, cnt in groups]
    for prod in product(*lists):
        flat = []
        for tup in prod:
            flat.extend(tup)
        yield tuple(flat)


def cotree_to_adj(tree):
    if tree == ("L",):
        return np.zeros((1, 1), dtype=int)
    op, children = tree
    mats = [cotree_to_adj(c) for c in children]
    sizes = [M.shape[0] for M in mats]
    n = sum(sizes)
    A = np.zeros((n, n), dtype=int)
    off = 0
    for M in mats:
        k = M.shape[0]
        A[off : off + k, off : off + k] = M
        off += k
    if op == "J":
        off_i = 0
        for si in sizes:
            off_j = 0
            for sj in sizes:
                if off_i != off_j:
                    A[off_i : off_i + si, off_j : off_j + sj] = 1
                off_j += sj
            off_i += si
    return A


def compute_nabla(incomp):
    k = incomp.shape[0]
    best = [0]
    def dfs(cand, sz):
        if sz + len(cand) <= best[0]:
            return
        if not cand:
            best[0] = max(best[0], sz)
            return
        v = cand[0]
        dfs([u for u in cand[1:] if incomp[v, u]], sz + 1)
        dfs(cand[1:], sz)
    dfs(list(range(k)), 0)
    return best[0]


def verify_table_2():
    print("[4/4] Verifying Table 2 (all 1,170 connected cographs, 2 <= n <= 9)...")
    connected = {1: [("L",)]}
    disconnected = {1: []}
    for n in range(2, 10):
        disconnected[n] = [("U", ch) for part in partitions(n, n - 1) for ch in gen_from_parts(connected, part)]
        co_conn = {k: ([("L",)] if k == 1 else disconnected[k]) for k in range(1, n)}
        connected[n] = [("J", ch) for part in partitions(n, n - 1) for ch in gen_from_parts(co_conn, part)]

    x = sp.Symbol("x")
    total_graphs = 0
    for n in range(2, 10):
        q_viols, dq_orig, dq_rep = 0, 0, 0
        for tree in connected[n]:
            total_graphs += 1
            A = cotree_to_adj(tree)
            N_open = [frozenset(np.where(A[u] == 1)[0]) for u in range(n)]
            N_closed = [N_open[u] | {u} for u in range(n)]
            visited = [False] * n
            classes = []
            for u in range(n):
                if visited[u]:
                    continue
                cls = [u]
                visited[u] = True
                for v in range(u + 1, n):
                    if not visited[v] and (N_open[u] - {v}) == (N_open[v] - {u}):
                        cls.append(v)
                        visited[v] = True
                classes.append(cls)
            k = len(classes)
            reps = [c[0] for c in classes]
            sizes = [len(c) for c in classes]
            is_clique = [bool(len(c) >= 2 and A[c[0], c[1]] == 1) for c in classes]
            deg = [int(np.sum(A[r])) for r in reps]
            comp = np.zeros((k, k), dtype=bool)
            for i in range(k):
                for j in range(k):
                    if N_open[reps[i]].issubset(N_closed[reps[j]]):
                        comp[i, j] = True
            incomp = ~(comp | comp.T)
            np.fill_diagonal(incomp, False)
            nabla = compute_nabla(incomp)
            T_Q = {deg[i] - (1 if is_clique[i] else 0) for i in range(k) if sizes[i] >= 2}
            T_DQ = {2 * (n - 1) - deg[i] - (1 if is_clique[i] else 2) for i in range(k) if sizes[i] >= 2}

            Qq = np.zeros((k, k), dtype=float)
            DQq = np.zeros((k, k), dtype=float)
            for i in range(k):
                Tr_i = 2 * (n - 1) - deg[i]
                Qq[i, i] = deg[i] + ((sizes[i] - 1) if is_clique[i] else 0)
                DQq[i, i] = Tr_i + (sizes[i] - 1) * (1 if is_clique[i] else 2)
                for j in range(k):
                    if i != j:
                        adj = int(A[reps[i], reps[j]])
                        Qq[i, j] = sizes[j] * adj
                        DQq[i, j] = sizes[j] * (1 if adj == 1 else 2)
            sq = np.sqrt(sizes)
            eQ = np.sort(np.linalg.eigvalsh((sq[:, None] * Qq) / sq[None, :]))
            eDQ = np.sort(np.linalg.eigvalsh((sq[:, None] * DQq) / sq[None, :]))
            if k >= 2 and np.min(np.diff(eQ)) < 1e-3:
                Mq = sp.Matrix([[int(round(Qq[i, j])) for j in range(k)] for i in range(k)])
                p = sp.Poly(Mq.charpoly(x).as_expr(), x)
                if sp.gcd(p, p.diff()).degree() > 0:
                    for e, mult in Mq.eigenvals().items():
                        if e not in T_Q and mult > nabla:
                            q_viols += 1
                            break
            if k >= 2 and np.min(np.diff(eDQ)) < 1e-3:
                Mdq = sp.Matrix([[int(round(DQq[i, j])) for j in range(k)] for i in range(k)])
                p = sp.Poly(Mdq.charpoly(x).as_expr(), x)
                if sp.gcd(p, p.diff()).degree() > 0:
                    for e, mult in Mdq.eigenvals().items():
                        if e not in T_DQ:
                            if mult > nabla:
                                dq_orig += 1
                            if mult > nabla + 1:
                                dq_rep += 1
        assert q_viols == 0 and dq_rep == 0
        print(f"      n={n}: {len(connected[n]):3d} graphs | Q>nabla: {q_viols} | DQ>nabla: {dq_orig} | DQ>nabla+1: {dq_rep}")
    assert total_graphs == 1170
    print(f"      PASSED: All {total_graphs} connected cographs verified.")


if __name__ == "__main__":
    verify_theorem_3_1()
    verify_proposition_3_2()
    verify_theorem_5_4()
    verify_table_2()
    print("\nALL VERIFICATIONS PASSED.")