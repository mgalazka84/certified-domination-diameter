"""Finite checks for certified domination and the smaller part.

Run without arguments to reproduce all results. Only the reduction checks
requires NumPy and SciPy. Vertices are zero-based. A graph is a sequence
of integers in which bit j of graph[i] records adjacency to j.
"""

import argparse
from collections import Counter, deque
from functools import cache
from itertools import combinations, combinations_with_replacement, permutations
from math import factorial, inf, prod


def indices(mask):
    """Yield the zero-based indices of the set bits."""
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def bipartite_graph(a, columns):
    """The first a vertices form A. Each column lists the neighbors in A."""
    graph = [0] * a + list(columns)
    for v, column in enumerate(columns, a):
        for x in indices(column):
            graph[x] |= 1 << v
    return graph


def distances(graph, start):
    result = [-1] * len(graph)
    result[start] = 0
    todo = deque([start])
    while todo:
        v = todo.popleft()
        for w in indices(graph[v]):
            if result[w] < 0:
                result[w] = result[v] + 1
                todo.append(w)
    return result


def connected_diameter(graph):
    diameter = 0
    for v in range(len(graph)):
        row = distances(graph, v)
        if -1 in row:
            return inf
        diameter = max(diameter, max(row))
    return diameter


def certified(graph, chosen):
    outside = ((1 << len(graph)) - 1) ^ chosen
    dominated = chosen
    for v in indices(chosen):
        if (graph[v] & outside).bit_count() == 1:
            return False
        dominated |= graph[v]
    return (dominated & outside) == outside


@cache
def choices_of_size(n, size):
    return tuple(sum(1 << v for v in chosen) for chosen in combinations(range(n), size))


def exhaustive_minimum(graph, limit=None):
    """Return infinity if no solution has size at most limit."""
    if limit is None:
        limit = len(graph)
    for k in range(limit + 1):
        if any(certified(graph, chosen) for chosen in choices_of_size(len(graph), k)):
            return k
    return inf


def classification(a, columns):
    if a == 1:
        return len(columns) >= 2
    if a == 2:
        return all(sum(bool(mask & (1 << v)) for mask in columns) >= 2 for v in range(2))
    if a == 3:
        return all(columns.count(mask) >= 2 for mask in (0b011, 0b101, 0b110))
    return False


def column_orderings(columns):
    """Number of labelings of B represented by these sorted columns."""
    return factorial(len(columns)) // prod(factorial(k) for k in Counter(columns).values())


def check_classification():
    all_graphs = 0
    all_equalities = 0
    four_vertex_side = 0
    print("Small diameter, fixed labels in both parts")
    print("A B  graphs  equality")
    for a in range(1, 5):
        for b in range(a, (6 if a == 3 else 5) + 1):
            eligible = equalities = 0
            for columns in combinations_with_replacement(range(1, 1 << a), b):
                graph = bipartite_graph(a, columns)
                if connected_diameter(graph) > 3:
                    continue
                val = exhaustive_minimum(graph, a)
                assert (val == a) == classification(a, columns), (a, columns, val)
                assert a < 4 or val < a, (a, columns, val)
                weight = column_orderings(columns)
                eligible += weight
                if val == a:
                    equalities += weight
            print(f"{a} {b}  {eligible}  {equalities}")
            all_graphs += eligible
            all_equalities += equalities
            if a == 4:
                four_vertex_side += eligible
    print(f"Total: {all_graphs} graphs, {all_equalities} equalities.")
    print(f"Cases with a = 4: {four_vertex_side}.")


def minimum_by_small_part(a, columns):
    """Section 4, with dictionaries for the reachable states in each layer."""
    best = inf
    for chosen_a in range(1 << a):
        states = {(0,) * a: 0}
        for column in columns:
            next_states = {}
            selected_neighbors = column & chosen_a
            other_neighbors = column & ~chosen_a
            for state, cost in states.items():
                if selected_neighbors:
                    new = list(state)
                    for x in indices(selected_neighbors):
                        new[x] = min(2, new[x] + 1)
                    key = tuple(new)
                    next_states[key] = min(next_states.get(key, inf), cost)
                if other_neighbors.bit_count() != 1:
                    new = list(state)
                    for x in indices(other_neighbors):
                        new[x] = 1
                    key = tuple(new)
                    next_states[key] = min(next_states.get(key, inf), cost + 1)
            states = next_states
        for state, cost in states.items():
            if all((value != 1 if chosen_a & (1 << x) else value == 1)
                   for x, value in enumerate(state)):
                best = min(best, chosen_a.bit_count() + cost)
    return best


def check_small_part_algorithm():
    total = 0
    for a, b in ((1, 3), (2, 3), (3, 3), (3, 4)):
        for columns in combinations_with_replacement(range(1 << a), b):
            graph = bipartite_graph(a, columns)
            val = minimum_by_small_part(a, columns)
            assert val == exhaustive_minimum(graph), (a, columns, val)
            total += column_orderings(columns)
    print(f"Small-part algorithm vs exhaustive minimum: {total} labeled graphs, including disconnected inputs.")


def base_graphs(n):
    """One representative per isomorphism class with no isolates (n <= 5)."""
    edges = list(combinations(range(n), 2))
    number = {edge: i for i, edge in enumerate(edges)}
    renamings = [
        [number[tuple(sorted((p[u], p[v])))] for u, v in edges]
        for p in permutations(range(n))
    ]
    seen = set()
    for mask in range(1, 1 << len(edges)):
        if mask in seen:
            continue
        present = list(indices(mask))
        for rename in renamings:
            seen.add(sum(1 << rename[i] for i in present))
        graph_edges = [edges[i] for i in present]
        if len({v for edge in graph_edges for v in edge}) == n:
            yield graph_edges


def vertex_cover_number(n, edges):
    for size in range(n + 1):
        for candidate in choices_of_size(n, size):
            if all(candidate & ((1 << u) | (1 << v)) for u, v in edges):
                return size
    raise AssertionError("The full vertex set must cover every edge.")


def reduction(n, edges, k):
    m = len(edges)
    t, s = 4*m - n, 4*m - 1 - k
    a = n + t + s
    columns = [(1 << a) - 1]  # r sees every vertex of A
    for u, v in edges:
        columns.extend([((1 << u) | (1 << v))] * 4)
    columns.extend(1 << w for w in range(n + t, a))
    return a, bipartite_graph(a, columns)


def optimum_by_integer_program(graph):
    """Generic 0-1 model of the definition; no forced vertices or cover data."""
    import numpy as np
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import lil_matrix

    n = len(graph)
    matrix = lil_matrix((3*n, 2*n), dtype=float)
    lower = np.full(3*n, -np.inf)
    upper = np.full(3*n, np.inf)
    for v, mask in enumerate(graph):
        ns = list(indices(mask))
        d = len(ns)
        matrix[3*v, [v] + ns] = 1
        lower[3*v] = 1
        for w in ns:
            matrix[3*v + 1, w] = 1
            matrix[3*v + 2, w] = 1
        matrix[3*v + 1, n+v] = 2
        upper[3*v + 1] = d
        matrix[3*v + 2, v] = -d
        matrix[3*v + 2, n+v] = d
        lower[3*v + 2] = 0
    answer = milp(c=np.r_[np.ones(n), np.zeros(n)], integrality=np.ones(2*n),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(matrix.tocsr(), lower, upper),
                  options={"mip_rel_gap": 0.0})
    if not answer.success:
        raise RuntimeError(answer.message)
    witness = sum(1 << v for v, value in enumerate(answer.x[:n]) if value > 0.5)
    assert certified(graph, witness)
    assert abs(witness.bit_count() - answer.fun) < 1e-5
    return witness.bit_count()


def check_reduction():
    originals = outputs = largest = 0
    for n in range(2, 6):
        for edges in base_graphs(n):
            tau = vertex_cover_number(n, edges)
            originals += 1
            for k in range(n + 1):
                a, graph = reduction(n, edges, k)
                val = optimum_by_integer_program(graph)
                assert val-a == tau-k, (n, edges, k, val, a, tau)
                assert len(graph) - a == a + 1
                assert max(distances(graph, a)) == 2
                assert connected_diameter(graph) == 4
                largest = max(largest, len(graph))
                outputs += 1
    print(f"Exact reduction: {originals} base graphs, {outputs} target instances, largest output order {largest}.")

    a, graph = reduction(2, [(0, 1)], 1)
    minima = [s for s in choices_of_size(len(graph), a) if certified(graph, s)]
    assert exhaustive_minimum(graph) == 6 and len(minima) == 2
    print(f"One-edge example: {len(graph)} vertices, optimum 6, {len(minima)} minimum sets.")


def fixed_diameter_reduction(n, edges, k, diameter):
    """Reattach existing leaf pairs along a path to obtain a prescribed diameter."""
    if diameter < 4:
        raise ValueError("The construction requires diameter at least four.")
    length = diameter - 3
    edges = list(edges)
    while 4 * len(edges) - 1 - k < length:
        edges.append((n, n + 1))
        n += 2
        k += 1

    hub, graph = reduction(n, edges, k)
    arm_count = 4 * len(edges) - 1 - k
    first_support = hub - arm_count
    for support in range(first_support + 1, first_support + length):
        previous = support - 1
        graph[hub] &= ~(1 << support)
        graph[support] &= ~(1 << hub)
        graph[previous] |= 1 << support
        graph[support] |= 1 << previous
    return hub, graph


def check_fixed_diameters():
    """Measure both part sizes after recoloring and solve the unrestricted 0-1 model."""
    instances = 0
    for n in range(2, 6):
        for edges in base_graphs(n):
            tau = vertex_cover_number(n, edges)
            for k in range(n + 1):
                for diameter in range(5, 11):
                    hub, graph = fixed_diameter_reduction(n, edges, k, diameter)
                    levels = distances(graph, hub)
                    assert -1 not in levels
                    assert all(abs(levels[v] - levels[w]) == 1
                               for v, neighbors in enumerate(graph) for w in indices(neighbors))
                    smaller = sum(level % 2 for level in levels)
                    assert len(graph) - smaller == smaller + 1
                    assert connected_diameter(graph) == diameter
                    val = optimum_by_integer_program(graph)
                    assert val - smaller == tau - k, (n, edges, k, diameter, val, smaller)
                    instances += 1
    print(f"Fixed diameters 5 through 10: {instances} instances; all diameters, part sizes and optimum gaps agree.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("check", choices=("all", "diameter", "algorithm", "reduction", "fixed-diameters"),
                        nargs="?", default="all")
    mode = parser.parse_args().check
    if mode in ("all", "diameter"):
        check_classification()
    if mode in ("all", "algorithm"):
        check_small_part_algorithm()
    if mode in ("all", "reduction"):
        check_reduction()
    if mode in ("all", "fixed-diameters"):
        check_fixed_diameters()
