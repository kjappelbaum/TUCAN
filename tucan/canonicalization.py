from tucan.graph_utils import sort_molecule_by_attribute, attribute_sequence
import networkx as nx
from itertools import pairwise
from operator import gt, lt, eq
from collections import deque


def partition_molecule_by_attribute(m, attribute):
    m_sorted = sort_molecule_by_attribute(m, attribute)
    attribute_sequences = [
        attribute_sequence(a, m_sorted, attribute) for a in sorted(m_sorted)
    ]
    updated_partitions = [0]
    for i, j in pairwise(attribute_sequences):
        current_partition = updated_partitions[-1]
        if i != j:
            current_partition += 1
        updated_partitions.append(current_partition)
    nx.set_node_attributes(
        m_sorted, dict(zip(sorted(m_sorted), updated_partitions)), "partition"
    )
    return m_sorted


def refine_partitions(m):
    current_partitions = list(
        nx.get_node_attributes(
            sort_molecule_by_attribute(m, "partition"), "partition"
        ).values()
    )
    m_refined = partition_molecule_by_attribute(m, "partition")
    refined_partitions = list(nx.get_node_attributes(m_refined, "partition").values())

    while current_partitions != refined_partitions:
        yield m_refined
        current_partitions = refined_partitions
        m_refined = partition_molecule_by_attribute(m_refined, "partition")
        refined_partitions = list(
            nx.get_node_attributes(m_refined, "partition").values()
        )


def assign_canonical_labels(m, traversal_priorities=[lt, gt, eq]):
    """Re-label nodes such that each node is connected to neighbors with the
    smallest possible labels.
    """
    partitions = m.nodes.data("partition")
    labels_by_partition = _labels_by_partition(m)
    final_labels = {}
    atom_queue = deque([0])
    nx.set_node_attributes(m, False, "explored")

    while atom_queue:
        a = atom_queue.popleft()
        if m.nodes[a]["explored"]:
            continue

        # assign smallest label available in this partition
        final_labels[a] = labels_by_partition[partitions[a]].pop()
        m.nodes[a]["explored"] = True

        neighbors = list(m.neighbors(a))
        for priority in traversal_priorities:
            atom_queue.extend(
                [n for n in neighbors if priority(partitions[a], partitions[n])]
            )

    nx.set_node_attributes(m, False, "explored")

    return final_labels


def _labels_by_partition(m):
    """Create dictionary of partitions to node labels."""
    partitions = set(sorted([v for _, v in m.nodes.data("partition")]))
    labels_by_partition = {p: set() for p in partitions}
    for a in m:
        labels_by_partition[m.nodes[a]["partition"]].add(a)
    labels_by_partition.update(
        (k, sorted(list(v), reverse=True)) for k, v in labels_by_partition.items()
    )
    return labels_by_partition


def _add_invariant_code(m):
    """Assign an invariant code to each atom (mutates graph)."""
    atomic_numbers = list(nx.get_node_attributes(m, "atomic_number").values())
    nx.set_node_attributes(
        m, dict(zip(list(m.nodes), atomic_numbers)), "invariant_code"
    )


def canonicalize_molecule(m):
    _add_invariant_code(m)
    m_partitioned_by_invariant_code = partition_molecule_by_attribute(
        m, "invariant_code"
    )
    m_refined = list(refine_partitions(m_partitioned_by_invariant_code))
    m_partitioned = m_refined[-1] if m_refined else m_partitioned_by_invariant_code
    canonical_labels = assign_canonical_labels(m_partitioned)
    return nx.relabel_nodes(m_partitioned, canonical_labels, copy=True)


# def bfs_molecule(m, root_idx):
#     """Breadth-first search over atoms.
#     Note that NetworkX provides the same algorithm in `dfs_edges()`.
#     This (re-)implementation allows for controlling the branching behavior
#     during the molecule traversal.
#     m: NetworkX graph.
#     root_idx: atom at which to start traversal.
#     """
#     m.nodes[root_idx]["explored"] = True
#     atom_queue = deque([root_idx])
#     while atom_queue:
#         a = atom_queue.popleft()
#         for n in m.neighbors(a):
#             if m.nodes[n]["explored"]:
#                 continue
#             yield (a, n)
#             m.nodes[n]["explored"] = True
#             atom_queue.append(n)

# def dfs_molecule(m, root_idx):
#     """Depth-first search over atoms.
#     Note that NetworkX provides the same algorithm in `bfs_edges()`.
#     This (re-)implementation allows for controlling the branching behavior
#     during the molecule traversal.
#     m: NetworkX graph.
#     root_idx: atom at which to start traversal.
#     """
#     m.nodes[root_idx]["explored"] = True
#     for n_idx in m.neighbors(root_idx):
#         if m.nodes[n_idx]["explored"]:
#             continue
#         yield (root_idx, n_idx)
#         yield from dfs_molecule(m, n_idx)

# def edge_dfs_molecule(m, root_idx):
#     """Depth-first search over edges.
#     Note that NetworkX provides the same algorithm in `edge_dfs ()`.
#     This (re-)implementation allows for controlling the branching behavior
#     during the molecule traversal.
#     m: NetworkX graph.
#     root_idx: atom at which to start traversal.
#     """
#     visited_edges = set()
#     visited_nodes = set()
#     edges = {}

#     nodes = list(m.nbunch_iter(root_idx))
#     for start_node in nodes:
#         stack = [start_node]
#         while stack:
#             current_node = stack[-1]
#             if current_node not in visited_nodes:
#                 edges[current_node] = iter(m.edges(current_node))
#                 visited_nodes.add(current_node)

#             try:
#                 edge = next(edges[current_node])
#             except StopIteration:
#                 # No more edges from the current node.
#                 stack.pop()
#             else:
#                 edgeid = (frozenset(edge[:2]),) + edge[2:]
#                 if edgeid not in visited_edges:
#                     visited_edges.add(edgeid)
#                     # Mark the traversed "to" node as to-be-explored.
#                     stack.append(edge[1])
#                     yield edge
