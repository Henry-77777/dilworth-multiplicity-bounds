#!/usr/bin/env python3
"""
100% Pure Exact Symbolic Verification Script (Zero Floating-Point Arithmetic)
accompanying the final manuscript:
  "A Spectral Dichotomy for Dilworth Multiplicity Bounds of the Signless Laplacian
   and Distance Signless Laplacian on C-Graphs and Cographs"

Scope & Reproducibility Note:
  All universal upper-bound theorems in the manuscript (Theorems 2.1, 4.4, 4.5, 5.3,
  and 5.8) are proved purely analytically; finite-order symbolic computations below
  confirm exact small-order enumerations and independently cross-check our explicit
  constructions, and do NOT substitute for analytical proofs of universal bounds.

What this script verifies using pure integer / rational / polynomial arithmetic (SymPy):
  1. Theorem 3.1: Infinite counterexample family G_k = C(k, 1, ..., 1) for k = 1..5
     on the full n x n integer matrix D^Q(G_k) using exact SymPy Matrix.rank().
  2. Proposition 3.2: Exhaustive pure symbolic irreducible-factor decomposition
     (sympy.factor_list over Q[x], zero root approximation) across all 127 even-length
     C-graph creation sequences of order 2 <= n <= 8, confirming H_8 = C(1^8) is the
     unique smallest counterexample to Mandal-Manna Conjecture 5.5 up to isomorphism.
  3. Theorem 5.9: Full n x n exact SymPy verification of B = K_36 V (3 * H_8)
     (n=60, nabla=3, exact 60x60 nullspace dimension of Q(B) - 40*I_60 equals 4) and
     F_s = K_1 V (s * B) for s = 2, 3 (exact integer eigenvector equations + exact SymPy rank).
  4. Table 2 (Remark 5.10): Exhaustive pure symbolic irreducible-factor decomposition
     (sympy.factor_list over Q[x], zero floating-point prefiltering, zero root solving)
     across ALL 1,170 non-isomorphic connected cographs of order 2 <= n <= 9, reporting
     exact counts of violating GRAPHS.

Usage:
  python3 verify_all.py
Requirements:
  Python 3.8+, numpy (integer array construction only), sympy
"""

from itertools import combinations_with_replacement, product
import platform
import random
import sys
import time
import numpy as np
import sympy as sp


def has_nonlocal_multiplicity_gt(poly, local_int_set, threshold):
    """
    Pure algebraic check over Q[x] via irreducible factorization (factor_list):
    Returns (True, witness_list) iff `poly` has at least one root not in `local_int_set`
    with algebraic multiplicity strictly greater than `threshold`.
    Because irreducible factors over Q are separable, each root of an irreducible factor
    `fac` has exact multiplicity equal to its exponent `exp`, with zero numerical root solving.
    """
    witnesses = []
    if sp.gcd(poly, poly.diff()).degree() == 0:
        # Square-free polynomial: every root has multiplicity 1
        if threshold >= 1:
            return False, witnesses
    _, factors = sp.factor_list(poly)
    for fac_poly, exp in factors:
        if exp <= threshold:
            continue
        fac = sp.Poly(fac_poly, poly.gen)
        if fac.degree() == 1:
            a, b = fac.all_coeffs()
            root = sp.Rational(-b, a)
            if root.is_integer and int(root) in local_int_set:
                continue
            witnesses.append((root, exp))
        else:
            # Irreducible of degree >= 2 over Q: all roots are irrational, hence not in local_int_set
            witnesses.append((fac.as_expr(), exp))
    return (len(witnesses) > 0), witnesses


def build_cgraph_recursive(alphas):
    """Build exact integer adjacency matrix of C(alpha_1, ..., alpha_m) from the recursive definition."""
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
    print("[1/5] Verifying Theorem 3.1 (infinite family G_k, k = 1..5, exact integer & SymPy rank)...")
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
        V_exact = sp.Matrix(np.column_stack([x1, x2]).tolist())
        assert V_exact.rank() == 2
    print("      PASSED: Exact integer eigenvectors and exact SymPy rank = 2 for k = 1..5.")


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
    print("[2/5] Verifying Proposition 3.2 (all 127 even-length C-graph sequences, 2 <= n <= 8, factor_list)...")
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
            viol, wits = has_nonlocal_multiplicity_gt(p, TQ, nabla)
            if viol:
                violations.append((alphas, nabla, wits))
    assert total == 127
    assert violations == [((1, 1, 1, 1, 1, 1, 1, 1), 1, [(sp.Rational(8), 2)])]
    print(f"      PASSED: {total} creation sequences checked via factor_list; unique violation is C(1^8) = H_8 (theta=8, mult=2 > nabla=1).")


def verify_theorem_5_9():
    print("[3/5] Verifying Theorem 5.9 (exact 60x60 nullspace of B & exact SymPy rank on F_2, F_3)...")
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
    assert sp.Matrix(np.column_stack(vecsB).tolist()).rank() == 4

    # Independent exact 60x60 SymPy nullspace confirming mult_{Q(B)}(40) == 4 EXACTLY
    QB_sp = sp.Matrix(QB.tolist())
    assert len((QB_sp - 40 * sp.eye(nB)).nullspace()) == 4

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
        assert sp.Matrix(np.column_stack(vecsF).tolist()).rank() == 4 * s - 1
    print("      PASSED: Exact 60x60 nullity(Q(B)-40I) = 4 and exact SymPy ranks on F_2 (7), F_3 (11).")


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
    print("[4/5] Verifying Table 2 / Remark 5.10 (ALL 1,170 connected cographs, 2 <= n <= 9, pure factor_list)...")
    t0 = time.time()
    connected = {1: [("L",)]}
    disconnected = {1: []}
    for n in range(2, 10):
        disconnected[n] = [("U", ch) for part in partitions(n, n - 1) for ch in gen_from_parts(connected, part)]
        co_conn = {k: ([("L",)] if k == 1 else disconnected[k]) for k in range(1, n)}
        connected[n] = [("J", ch) for part in partitions(n, n - 1) for ch in gen_from_parts(co_conn, part)]

    x = sp.Symbol("x")
    expected_graphs = {2: 1, 3: 2, 4: 5, 5: 12, 6: 33, 7: 90, 8: 261, 9: 766}
    total_graphs = 0
    for n in range(2, 10):
        q_viol_graphs, dq_orig_graphs, dq_rep_graphs = 0, 0, 0
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

            Mq = sp.zeros(k, k)
            Mdq = sp.zeros(k, k)
            for i in range(k):
                Tr_i = 2 * (n - 1) - deg[i]
                Mq[i, i] = deg[i] + ((sizes[i] - 1) if is_clique[i] else 0)
                Mdq[i, i] = Tr_i + (sizes[i] - 1) * (1 if is_clique[i] else 2)
                for j in range(k):
                    if i != j:
                        adj = int(A[reps[i], reps[j]])
                        Mq[i, j] = sizes[j] * adj
                        Mdq[i, j] = sizes[j] * (1 if adj == 1 else 2)

            pQ = sp.Poly(Mq.charpoly(x).as_expr(), x)
            if has_nonlocal_multiplicity_gt(pQ, T_Q, nabla)[0]:
                q_viol_graphs += 1

            pDQ = sp.Poly(Mdq.charpoly(x).as_expr(), x)
            violates, witnesses = has_nonlocal_multiplicity_gt(pDQ, T_DQ, nabla)
            if violates:
                # One rational witness represents one eigenvalue; a nonlinear
                # irreducible factor would represent multiple distinct roots.
                assert len(witnesses) == 1 and witnesses[0][0].is_Rational, (tree, witnesses)
                dq_orig_graphs += 1
            if has_nonlocal_multiplicity_gt(pDQ, T_DQ, nabla + 1)[0]:
                dq_rep_graphs += 1

        observed = (len(connected[n]), q_viol_graphs, dq_orig_graphs, dq_rep_graphs)
        expected = (expected_graphs[n], 0, 4 if n == 8 else 0, 0)
        assert observed == expected, f"Table 2, n={n}: expected {expected}, got {observed}"
        print(
            f"      n={n}: {len(connected[n]):3d} graphs | "
            f"Graphs(Q>nabla): {q_viol_graphs} | Graphs(DQ>nabla): {dq_orig_graphs} | Graphs(DQ>nabla+1): {dq_rep_graphs}"
        )
    assert total_graphs == 1170
    print(f"      PASSED: All {total_graphs} connected cographs verified via factor_list in {time.time()-t0:.2f}s.")




def verify_theorem_5_3_case_IV():
    print("[5/5] Verifying Theorem 5.3 Case (IV) (102 weighted threshold graphs via pure factor_list)...")
    x = sp.Symbol("x")

    def check_case_IV(b_sizes, a_sizes, wB, wA):
        s = len(b_sizes)
        assert s >= 2 and b_sizes[0] == 1 and len(a_sizes) == s - 1
        blocks = [("B", 1, b_sizes[0], wB[0])]
        for k in range(2, s + 1):
            blocks.append(("A", k, a_sizes[k - 2], wA[k - 2]))
            blocks.append(("B", k, b_sizes[k - 1], wB[k - 1]))
        verts = []
        for typ, idx, sz, wt in blocks:
            for _ in range(sz):
                verts.append((typ, idx, wt))
        n = len(verts)
        A = sp.zeros(n, n)
        for u in range(n):
            tu, iu, _ = verts[u]
            for v in range(u + 1, n):
                tv, iv, _ = verts[v]
                adj = (
                    (tu == "B" and tv == "B")
                    or (tu == "A" and tv == "B" and iu <= iv)
                    or (tu == "B" and tv == "A" and iv <= iu)
                )
                if adj:
                    A[u, v] = 1
                    A[v, u] = 1
        deg = [sum(A.row(u)) for u in range(n)]
        Qw = sp.diag(*[deg[u] + verts[u][2] for u in range(n)]) + A
        T_local = set()
        off = 0
        for typ, idx, sz, wt in blocks:
            if sz >= 2:
                d = deg[off]
                T_local.add(int(d - 1 + wt) if typ == "B" else int(d + wt))
            off += sz
        p = sp.Poly(Qw.charpoly(x).as_expr(), x)
        viol, wits = has_nonlocal_multiplicity_gt(p, T_local, 1)
        assert not viol, f"Unexpected Case (IV) violation: {wits}"

    # Explicit star boundary examples: K_1 V (K_1 U K_1), a1=0, b1=1, a2=b2=1
    check_case_IV([1, 1], [1], [1, 2], [0])
    check_case_IV([1, 1], [1], [2, 2], [0])

    # 100 deterministic random Case (IV) instances (40 with wB1 - wA2 == 1, 60 with jump >= 1)
    rng = random.Random(42)
    for trial in range(100):
        s = rng.choice([2, 3, 4])
        b_sizes = [1] + [rng.randint(1, 2) for _ in range(s - 1)]
        a_sizes = [rng.randint(1, 2) for _ in range(s - 1)]
        wA_rev = [rng.randint(-2, 2)]
        for _ in range(s - 2):
            wA_rev.append(wA_rev[-1] + rng.randint(0, 2))
        wA = list(reversed(wA_rev))
        jump = 1 if trial < 40 else rng.randint(1, 4)
        wB = [wA[0] + jump]
        for _ in range(s - 1):
            wB.append(wB[-1] + rng.randint(0, 2))
        check_case_IV(b_sizes, a_sizes, wB, wA)
    print("      PASSED: 102/102 Case (IV) weighted threshold graphs have zero non-local repeated roots.")


if __name__ == "__main__":
    if not __debug__:
        sys.exit("Run without -O or -OO to enable verification assertions.")
    print(f"Environment: Python {platform.python_version()} | NumPy {np.__version__} | SymPy {sp.__version__}")
    print("-" * 78)
    verify_theorem_3_1()
    verify_proposition_3_2()
    verify_theorem_5_9()
    verify_table_2()
    verify_theorem_5_3_case_IV()
    print("\nALL 100% PURE SYMBOLIC VERIFICATIONS PASSED.")