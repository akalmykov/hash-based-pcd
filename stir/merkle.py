from hashlib import sha3_256
from math import log2, ceil
from utils import is_power_of_two
import galois
from dataclasses import dataclass
from typing import List


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


class MerkleTree:

    def root(self):
        return self.non_leaf_nodes[0]

    def height(self):
        return int(log2(len(self.leaf_nodes))) + 1

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


# class MerkleTree:
#     """
#     A simple and naive implementation of an immutable Merkle tree.
#     """

#     def __init__(self, GF, data):
#         assert isinstance(data, list)
#         assert len(data) > 0, "Cannot construct an empty Merkle Tree."
#         num_leaves = 2 ** ceil(log2(len(data)))
#         self.data = data + [GF(0)] * (num_leaves - len(data))
#         self.height = int(log2(num_leaves))
#         self.facts = {}
#         self.root = self.build_tree()

#     def get_authentication_path(self, leaf_id):
#         assert 0 <= leaf_id < len(self.data)
#         node_id = leaf_id + len(self.data)
#         cur = self.root
#         decommitment = []
#         # In a Merkle Tree, the path from the root to a leaf, corresponds to the the leaf id's
#         # binary representation, starting from the second-MSB, where '0' means 'left', and '1' means
#         # 'right'.
#         # We therefore iterate over the bits of the binary representation - skipping the '0b'
#         # prefix, as well as the MSB.
#         for bit in bin(node_id)[3:]:
#             cur, auth = self.facts[cur]
#             if bit == "1":
#                 auth, cur = cur, auth
#             decommitment.append(auth)
#         return decommitment

#     def build_tree(self):
#         return self.recursive_build_tree(1)

#     def recursive_build_tree(self, node_id):
#         if node_id >= len(self.data):
#             # A leaf.
#             id_in_data = node_id - len(self.data)
#             leaf_data = str(self.data[id_in_data])
#             h = sha256(leaf_data.encode()).hexdigest()
#             self.facts[h] = leaf_data
#             return h
#         else:
#             # An internal node.
#             left = self.recursive_build_tree(node_id * 2)
#             right = self.recursive_build_tree(node_id * 2 + 1)
#             h = sha256((left + right).encode()).hexdigest()
#             self.facts[h] = (left, right)
#             return h


# def verify_decommitment(leaf_id, leaf_data, decommitment, root):
#     leaf_num = 2 ** len(decommitment)
#     node_id = leaf_id + leaf_num
#     cur = sha256(str(leaf_data).encode()).hexdigest()
#     for bit, auth in zip(bin(node_id)[3:][::-1], decommitment[::-1]):
#         if bit == "0":
#             h = cur + auth
#         else:
#             h = auth + cur
#         cur = sha256(h.encode()).hexdigest()
#     return cur == root
