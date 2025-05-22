from hashlib import sha3_256
from math import log2, ceil
from .utils import is_power_of_two
import galois
import pickle
from dataclasses import dataclass
from typing import List, TYPE_CHECKING


def serialize_field_vector(field_vector):
    vector_len = len(field_vector)
    serialized_bytes = vector_len.to_bytes(8, "little") + b"".join(
        [int(x).to_bytes(24, "little") for x in field_vector]
    )
    return serialized_bytes


def hash_field_vector(field_vector):
    sha3 = sha3_256()
    sha3.update(serialize_field_vector(field_vector))
    return list(sha3.digest())


def hash_pair(left_hash, right_hash):
    hasher = sha3_256()
    hasher.update(bytes(left_hash))
    hasher.update(bytes(right_hash))
    return list(hasher.digest())


@dataclass
class MerkleMultiProof:
    indexes: List[int]
    auth_paths_prefix_lengths: List[int]
    auth_paths_suffixes: List[List[bytes]]
    leaf_siblings_hashes: List[bytes]

    def verify(self, root_hash, leaves):
        """Verify that leaves at self.indexes of the merkle tree match the given root hash."""
        if len(self.indexes) == 0:
            return True

        def prefix_decode_path(prev_path, prefix_length, suffix):
            """Decode a path using prefix encoding."""
            return prev_path[:prefix_length] + suffix

        def convert_index_to_last_level(index, tree_height):
            """Convert index to the last level index."""
            return index + (1 << (tree_height - 1)) - 1

        def parent(index):
            """Get the parent index."""
            return (index - 1) >> 1

        def select_left_right_child(index, curr_hash, sibling_hash):
            """Select left and right child based on index."""
            if index % 2 == 0:
                return curr_hash, sibling_hash
            else:
                return sibling_hash, curr_hash

        tree_height = len(self.auth_paths_suffixes[0]) + 2

        # Lookup table to avoid redundant hash computations
        hash_lut = {}

        # Init prev path for decoding
        prev_path = self.auth_paths_suffixes[0]

        for i in range(len(self.indexes)):
            leaf_index = self.indexes[i]
            leaf = leaves[i]
            leaf_sibling_hash = self.leaf_siblings_hashes[i]

            # Decode i-th auth path
            auth_path = prefix_decode_path(
                prev_path,
                self.auth_paths_prefix_lengths[i],
                self.auth_paths_suffixes[i],
            )
            # Update prev path for decoding next one
            prev_path = auth_path

            claimed_leaf_hash = hash_field_vector(leaf)

            # Determine if the leaf is the left or right child
            left_child, right_child = select_left_right_child(
                leaf_index, claimed_leaf_hash, leaf_sibling_hash
            )

            # Track position in the path
            index = leaf_index
            index_in_tree = convert_index_to_last_level(leaf_index, tree_height)
            index >>= 1
            index_in_tree = parent(index_in_tree)

            # Compute hash of the leaf level
            curr_path_node = hash_lut.get(index_in_tree)
            if curr_path_node is None:
                curr_path_node = hash_pair(left_child, right_child)
                hash_lut[index_in_tree] = curr_path_node

            # Check levels between leaf level and root
            for level in range(len(auth_path) - 1, -1, -1):
                # Determine if current node is left or right child
                left, right = select_left_right_child(
                    index, curr_path_node, auth_path[level]
                )

                # Update position
                index >>= 1
                index_in_tree = parent(index_in_tree)

                # Compute hash at this level
                curr_path_node = hash_lut.get(index_in_tree)
                if curr_path_node is None:
                    curr_path_node = hash_pair(left, right)
                    hash_lut[index_in_tree] = curr_path_node

            # Check if final hash is root
            if bytes(curr_path_node) != bytes(root_hash):
                return False

        return True


class MerkleTree:

    def root(self):
        return self.non_leaf_nodes[0]

    def height(self):
        return int(log2(len(self.leaf_nodes))) + 1

    def serialize(self, filepath: str) -> None:
        with open(filepath, "wb") as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(filepath: str) -> "MerkleTree":
        with open(filepath, "rb") as f:
            tree = pickle.load(f)
        return tree

    def __init__(self, leaf_nodes_data):
        size_leaf = len(leaf_nodes_data)
        assert size_leaf > 1 and is_power_of_two(size_leaf)
        self.leaf_nodes = [hash_field_vector(x) for x in leaf_nodes_data]
        # Calculate non-leaf nodes from stacked_evals
        # Store the non-leaf nodes in level order. The first element is the root node.
        # The ith node's (starting at index 0) children are at indices 2*i+1, 2*i+2
        size_non_leaf = size_leaf - 1
        # hash_of_empty = int(0).to_bytes(32)
        self.non_leaf_nodes = [None] * size_non_leaf

        for i, j in zip(
            range(size_leaf - 1, 0, -2),
            range(size_non_leaf - 1, 0, -1),
        ):
            self.non_leaf_nodes[j] = hash_pair(
                self.leaf_nodes[i - 1], self.leaf_nodes[i]
            )

        for i in range(j - 1, -1, -1):
            left_child_index = 2 * i + 1
            right_child_index = 2 * i + 2
            left_child = self.non_leaf_nodes[left_child_index]
            right_child = self.non_leaf_nodes[right_child_index]
            self.non_leaf_nodes[i] = hash_pair(left_child, right_child)

    def get_leaf_sibling_hash(self, index):
        if index & 1 == 0:
            return self.leaf_nodes[index + 1]
        else:
            return self.leaf_nodes[index - 1]

    def convert_index_to_last_level(self, index, tree_height):
        return index + (1 << (tree_height - 1)) - 1

    def parent(self, index):
        return (index - 1) >> 1

    def is_left_non_leaf_child(self, index):
        return index & 1 == 1

    def non_leaf_sibling(self, index):
        if index == 0:
            return None
        elif self.is_left_non_leaf_child(index):
            return index + 1
        else:
            return index - 1

    def compute_auth_path(self, index):
        tree_height = self.height()
        leaf_index_in_tree = self.convert_index_to_last_level(index, tree_height)
        # len(path) = `tree height - 2`, the two missing elements being the leaf sibling hash and the root
        path = []
        current_node = self.parent(leaf_index_in_tree)
        while current_node != 0:
            sibling_node = self.non_leaf_sibling(current_node)
            path.append(self.non_leaf_nodes[sibling_node])
            current_node = self.parent(current_node)

        assert len(path) == tree_height - 2
        path.reverse()
        return path

    def prefix_encode_path(self, prev_path, path):
        prefix_length = 0
        for a, b in zip(prev_path, path):
            if a != b:
                break
            prefix_length += 1

        return (prefix_length, path[prefix_length:])

    def generate_multi_proof(self, indexes):
        indexes = sorted(list(set(indexes)))
        auth_paths_prefix_lengths = []
        auth_paths_suffixes = []
        leaf_siblings_hashes = []
        prev_path = []
        for index in indexes:
            leaf_siblings_hashes.append(self.get_leaf_sibling_hash(index))
            path = self.compute_auth_path(index)
            prefix_len, suffix = self.prefix_encode_path(prev_path, path)
            auth_paths_prefix_lengths.append(prefix_len)
            auth_paths_suffixes.append(suffix)
            prev_path = path

        return MerkleMultiProof(
            indexes=indexes,
            auth_paths_prefix_lengths=auth_paths_prefix_lengths,
            auth_paths_suffixes=auth_paths_suffixes,
            leaf_siblings_hashes=leaf_siblings_hashes,
        )
