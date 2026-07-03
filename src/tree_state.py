"""
Phase 1: TreeState dataclass representing T = (V, E, r, l, X).
Tree-valued state space for evolutionary processes.
"""
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class TreeState:
    """
    Represents a phylogenetic tree state T = (V, E, r, l, X).

    Attributes:
        node_ids: ordered list of node identifiers (V)
        root_id: root node identifier (r)
        edges: list of (parent, child) tuples (E)
        branch_lengths: dict mapping (parent, child) → length (l)
        node_seqs: dict mapping node_id → sequence string (X)
        active_leaves: list of leaf node IDs that are currently growing
    """
    node_ids: list[str]
    root_id: str
    edges: list[tuple[str, str]]
    branch_lengths: dict[tuple[str, str], float]
    node_seqs: dict[str, str]
    active_leaves: list[str] = field(default_factory=list)

    def __post_init__(self):
        """Validate tree consistency."""
        if self.root_id not in self.node_ids:
            raise ValueError(f"Root {self.root_id} not in node_ids")
        if not all(seq_id in self.node_ids for seq_id in self.node_seqs.keys()):
            raise ValueError("node_seqs contains unknown nodes")
        if not all(p in self.node_ids and c in self.node_ids for p, c in self.edges):
            raise ValueError("edges reference unknown nodes")

    def n_nodes(self) -> int:
        """Total number of nodes."""
        return len(self.node_ids)

    def n_leaves(self) -> int:
        """Number of leaf nodes (nodes with no children)."""
        return sum(1 for node_id in self.node_ids if self.is_leaf(node_id))

    def get_children(self, node_id: str) -> list[str]:
        """Get all direct children of a node."""
        return [c for p, c in self.edges if p == node_id]

    def get_parent(self, node_id: str) -> Optional[str]:
        """Get parent of a node, or None if root."""
        for p, c in self.edges:
            if c == node_id:
                return p
        return None

    def get_siblings(self, node_id: str) -> list[str]:
        """Get sibling nodes (same parent)."""
        parent = self.get_parent(node_id)
        if parent is None:
            return []
        return [c for c in self.get_children(parent) if c != node_id]

    def is_leaf(self, node_id: str) -> bool:
        """Check if node is a leaf."""
        return len(self.get_children(node_id)) == 0

    def apply_mutation(self, node_id: str, position: int, aa: str) -> "TreeState":
        """Apply single amino acid mutation at position in node sequence."""
        if node_id not in self.node_seqs:
            raise ValueError(f"Node {node_id} not found")

        seq = self.node_seqs[node_id]
        if position < 0 or position >= len(seq):
            raise ValueError(f"Position {position} out of range [0, {len(seq)})")

        # Create new sequence with mutation
        new_seq = seq[:position] + aa + seq[position + 1 :]
        new_node_seqs = self.node_seqs.copy()
        new_node_seqs[node_id] = new_seq

        return TreeState(
            node_ids=self.node_ids.copy(),
            root_id=self.root_id,
            edges=self.edges.copy(),
            branch_lengths=self.branch_lengths.copy(),
            node_seqs=new_node_seqs,
            active_leaves=self.active_leaves.copy(),
        )

    def branch_node(self, node_id: str, child_seqs: list[str]) -> "TreeState":
        """
        Create branching event: node produces multiple children with given sequences.
        Children are named as f"{node_id}_child_{i}".
        """
        if node_id not in self.node_seqs:
            raise ValueError(f"Node {node_id} not found")

        if len(child_seqs) == 0:
            raise ValueError("Must create at least one child")

        # Validate sequence lengths match
        seq_len = len(self.node_seqs[node_id])
        if not all(len(seq) == seq_len for seq in child_seqs):
            raise ValueError("All child sequences must have same length as parent")

        new_node_ids = self.node_ids.copy()
        new_edges = self.edges.copy()
        new_node_seqs = self.node_seqs.copy()
        new_branch_lengths = self.branch_lengths.copy()

        for i, child_seq in enumerate(child_seqs):
            child_id = f"{node_id}_child_{i}"
            new_node_ids.append(child_id)
            new_edges.append((node_id, child_id))
            new_node_seqs[child_id] = child_seq
            new_branch_lengths[(node_id, child_id)] = 0.0

        # Remove parent from active leaves, add children
        new_active_leaves = [n for n in self.active_leaves if n != node_id]
        new_active_leaves.extend([f"{node_id}_child_{i}" for i in range(len(child_seqs))])

        return TreeState(
            node_ids=new_node_ids,
            root_id=self.root_id,
            edges=new_edges,
            branch_lengths=new_branch_lengths,
            node_seqs=new_node_seqs,
            active_leaves=new_active_leaves,
        )

    def extend_branch(self, node_id: str, delta: float) -> "TreeState":
        """Extend branch leading to node by delta time units."""
        if node_id == self.root_id:
            raise ValueError("Cannot extend branch leading to root")

        parent = self.get_parent(node_id)
        if parent is None:
            raise ValueError(f"Node {node_id} has no parent")

        key = (parent, node_id)
        if key not in self.branch_lengths:
            raise ValueError(f"No branch length for {key}")

        new_branch_lengths = self.branch_lengths.copy()
        new_branch_lengths[key] += delta

        return TreeState(
            node_ids=self.node_ids.copy(),
            root_id=self.root_id,
            edges=self.edges.copy(),
            branch_lengths=new_branch_lengths,
            node_seqs=self.node_seqs.copy(),
            active_leaves=self.active_leaves.copy(),
        )

    def terminate_leaf(self, node_id: str) -> "TreeState":
        """Mark a leaf as terminal (stop growing)."""
        if node_id not in self.active_leaves:
            raise ValueError(f"Node {node_id} is not an active leaf")

        new_active_leaves = [n for n in self.active_leaves if n != node_id]

        return TreeState(
            node_ids=self.node_ids.copy(),
            root_id=self.root_id,
            edges=self.edges.copy(),
            branch_lengths=self.branch_lengths.copy(),
            node_seqs=self.node_seqs.copy(),
            active_leaves=new_active_leaves,
        )

    def to_dict(self) -> dict:
        """Serialize to dictionary."""
        return {
            "node_ids": self.node_ids,
            "root_id": self.root_id,
            "edges": self.edges,
            "branch_lengths": {str(k): v for k, v in self.branch_lengths.items()},
            "node_seqs": self.node_seqs,
            "active_leaves": self.active_leaves,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TreeState":
        """Deserialize from dictionary."""
        # Convert string keys back to tuples for branch_lengths
        branch_lengths = {}
        for k, v in d["branch_lengths"].items():
            if isinstance(k, str):
                # Parse "(parent, child)" string format safely
                parts = k.strip("()").split(", ")
                if len(parts) == 2:
                    k = (parts[0].strip().strip("'\""), parts[1].strip().strip("'\""))
                else:
                    raise ValueError(f"Cannot parse branch_lengths key: {k}")
            branch_lengths[k] = v

        return cls(
            node_ids=d["node_ids"],
            root_id=d["root_id"],
            edges=d["edges"],
            branch_lengths=branch_lengths,
            node_seqs=d["node_seqs"],
            active_leaves=d.get("active_leaves", []),
        )

    @classmethod
    def root_only(cls, x0: str) -> "TreeState":
        """Create initial tree state with only root node."""
        root_id = "root"
        return cls(
            node_ids=[root_id],
            root_id=root_id,
            edges=[],
            branch_lengths={},
            node_seqs={root_id: x0},
            active_leaves=[root_id],
        )

    @classmethod
    def from_newick_pkl(cls, pkl_data: dict) -> "TreeState":
        """
        Create TreeState from preprocessed tree pickle (Phase 0 output).

        Expected pkl_data keys: name, root_id, root_seq, node_seqs, edges, branch_lengths, n_leaves, n_nodes
        """
        parent_ids = {parent for parent, _ in pkl_data["edges"]}
        active_leaves = [
            node_id
            for node_id in pkl_data["node_seqs"]
            if node_id not in parent_ids
        ]
        return cls(
            node_ids=list(pkl_data["node_seqs"].keys()),
            root_id=pkl_data["root_id"],
            edges=pkl_data["edges"],
            branch_lengths={tuple(e): pkl_data["branch_lengths"][e] for e in pkl_data["edges"]},
            node_seqs=pkl_data["node_seqs"],
            active_leaves=active_leaves,
        )
