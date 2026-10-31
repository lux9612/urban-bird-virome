                                                                              

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass(eq=False)
class Node:
    name: str = ""
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.position = 0

    def parse(self) -> Node:
        root = self._subtree()
        self._skip_space()
        if self._peek() == ";":
            self.position += 1
        self._skip_space()
        if self.position != len(self.text):
            raise ValueError(f"Unexpected Newick text at position {self.position}")
        root.parent = None
        root.length = 0.0
        return root

    def _subtree(self) -> Node:
        self._skip_space()
        if self._peek() == "(":
            self.position += 1
            children = [self._subtree()]
            while True:
                self._skip_space()
                char = self._peek()
                if char == ",":
                    self.position += 1
                    children.append(self._subtree())
                elif char == ")":
                    self.position += 1
                    break
                else:
                    raise ValueError(f"Expected ',' or ')' at position {self.position}")
            name = self._label()
            length = self._length()
            node = Node(name=name, length=length, children=children)
            for child in children:
                child.parent = node
            return node
        name = self._label(required=True)
        return Node(name=name, length=self._length())

    def _label(self, required: bool = False) -> str:
        self._skip_space()
        if self._peek() == "'":
            self.position += 1
            start = self.position
            while self.position < len(self.text) and self.text[self.position] != "'":
                self.position += 1
            if self.position >= len(self.text):
                raise ValueError("Unterminated quoted Newick label")
            label = self.text[start : self.position]
            self.position += 1
            return label
        start = self.position
        while self.position < len(self.text) and self.text[self.position] not in ":,();":
            if not self.text[self.position].isspace():
                self.position += 1
            else:
                break
        label = self.text[start : self.position].strip()
        if required and not label:
            raise ValueError(f"Missing Newick tip label at position {start}")
        return label

    def _length(self) -> float:
        self._skip_space()
        if self._peek() != ":":
            return 0.0
        self.position += 1
        start = self.position
        while self.position < len(self.text) and self.text[self.position] not in ",();":
            self.position += 1
        token = self.text[start : self.position].strip()
        if not token:
            raise ValueError(f"Missing branch length at position {start}")
        return float(token)

    def _skip_space(self) -> None:
        while self.position < len(self.text) and self.text[self.position].isspace():
            self.position += 1

    def _peek(self) -> str:
        return self.text[self.position] if self.position < len(self.text) else ""


def read_newick(path: Path) -> Node:
    return NewickParser(path.read_text(encoding="utf-8-sig")).parse()


def leaves(root: Node) -> list[Node]:
    output: list[Node] = []
    stack = [root]
    while stack:
        node = stack.pop()
        if node.is_leaf:
            output.append(node)
        else:
            stack.extend(reversed(node.children))
    return output


def leaf_lookup(root: Node) -> dict[str, Node]:
    result = {leaf.name: leaf for leaf in leaves(root)}
    if len(result) != len(leaves(root)):
        raise ValueError("Newick leaf names are not unique")
    return result


def root_distance(node: Node) -> float:
    distance = 0.0
    current: Node | None = node
    while current is not None and current.parent is not None:
        distance += current.length
        current = current.parent
    return distance


def patristic_distance(first: Node, second: Node) -> float:
    first_ancestors: dict[Node, float] = {}
    distance = 0.0
    current: Node | None = first
    while current is not None:
        first_ancestors[current] = distance
        if current.parent is None:
            break
        distance += current.length
        current = current.parent

    distance = 0.0
    current = second
    while current not in first_ancestors:
        if current is None or current.parent is None:
            raise ValueError("Nodes do not share a common root")
        distance += current.length
        current = current.parent
    return first_ancestors[current] + distance


def distance_matrix(root: Node, ordered_names: Iterable[str]) -> list[list[float]]:
    lookup = leaf_lookup(root)
    names = list(ordered_names)
    return [[patristic_distance(lookup[a], lookup[b]) for b in names] for a in names]


def prune_and_rename(root: Node, selected_to_output: dict[str, str]) -> Node:
                                                                                    

    def clone(node: Node) -> Node | None:
        if node.is_leaf:
            if node.name not in selected_to_output:
                return None
            return Node(name=selected_to_output[node.name], length=node.length)
        retained = [child for child in (clone(child) for child in node.children) if child is not None]
        if not retained:
            return None
        if len(retained) == 1:
            child = retained[0]
            child.length += node.length
            return child
        copied = Node(name="", length=node.length, children=retained)
        for child in retained:
            child.parent = copied
        return copied

    pruned = clone(root)
    if pruned is None:
        raise ValueError("Pruning removed every tree tip")
    pruned.parent = None
    pruned.length = 0.0
    return pruned


def to_newick(root: Node) -> str:
    def format_node(node: Node, is_root: bool = False) -> str:
        if node.is_leaf:
            label = node.name
        else:
            label = "(" + ",".join(format_node(child) for child in node.children) + ")"
        if not is_root:
            label += f":{node.length:.8f}"
        return label

    return format_node(root, is_root=True) + ";\n"

