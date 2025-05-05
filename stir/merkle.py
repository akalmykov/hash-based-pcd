from hashlib import sha3_256
from math import log2, ceil
from utils import is_power_of_two
import galois


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


class MerkleTree:

    def height(self):
        ceil(log2(len(self.leaf_nodes)))

    def __init__(self, leaf_nodes_data):
        size_leaf = len(leaf_nodes_data)
        assert size_leaf > 1 and is_power_of_two(size_leaf)
        self.leaf_nodes = [hash_field_vector(x) for x in leaf_nodes_data]
        # Calculate non-leaf nodes from stacked_evals
        # Store the non-leaf nodes in level order. The first element is the root node.
        # The ith node's (starting at index 0) children are at indices 2*i+1, 2*i+2
        tree_height = self.height()
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
