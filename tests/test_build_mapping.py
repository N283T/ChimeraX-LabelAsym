"""Unit tests for ``_build_mapping`` in ``hook.py``.

These tests exercise the pure-Python mapping logic with a fake CIFTable so
they do not require a ChimeraX runtime.
"""

import hook


class FakeAtomSite:
    """Minimal stand-in for ``chimerax.mmcif.CIFTable``."""

    def __init__(self, columns, rows):
        self._columns = list(columns)
        self._rows = rows

    def has_field(self, name):
        return name in self._columns

    def fields(self, requested):
        idx = [self._columns.index(f) for f in requested]
        return [[row[i] for i in idx] for row in self._rows]


def test_basic_auth_keyed_mapping():
    atom_site = FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id"],
        rows=[("E", "A", "142")],
    )
    mapping, stats = hook._build_mapping(atom_site)

    assert mapping == {("auth", "A", 142, ""): "E"}
    assert stats["skipped_auth_seq"] == 0
    assert stats["missing_fields"] == []


def test_insertion_code_normalization():
    atom_site = FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"],
        rows=[
            ("A", "A", "1", "?"),
            ("A", "A", "2", "."),
            ("A", "A", "3", ""),
            ("A", "A", "4", "B"),
        ],
    )
    mapping, _ = hook._build_mapping(atom_site)

    assert ("auth", "A", 1, "") in mapping
    assert ("auth", "A", 2, "") in mapping
    assert ("auth", "A", 3, "") in mapping
    assert ("auth", "A", 4, "B") in mapping


def test_duplicate_rows_first_wins():
    """atom_site has many rows per residue; setdefault must keep the first."""
    atom_site = FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id"],
        rows=[("E", "A", "1"), ("X", "A", "1"), ("Y", "A", "1")],
    )
    mapping, _ = hook._build_mapping(atom_site)

    assert mapping == {("auth", "A", 1, ""): "E"}


def test_non_integer_auth_seq_counted():
    atom_site = FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id"],
        rows=[
            ("A", "A", "1"),
            ("E", "A", "."),
            ("E", "A", "?"),
        ],
    )
    mapping, stats = hook._build_mapping(atom_site)

    assert mapping == {("auth", "A", 1, ""): "A"}
    assert stats["skipped_auth_seq"] == 2


def test_label_seq_id_adds_label_keyed_entry():
    atom_site = FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id", "label_seq_id"],
        rows=[("A", "A", "10", "1")],
    )
    mapping, stats = hook._build_mapping(atom_site)

    assert ("auth", "A", 10, "") in mapping
    assert ("label", "A", 1, "") in mapping
    assert stats["skipped_label_seq"] == 0


def test_label_seq_id_skipped_when_dot():
    atom_site = FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id", "label_seq_id"],
        rows=[("E", "A", "142", ".")],
    )
    mapping, stats = hook._build_mapping(atom_site)

    assert ("auth", "A", 142, "") in mapping
    assert ("label", "E", 0, "") not in mapping
    assert stats["skipped_label_seq"] == 1


def test_missing_required_field_reports_and_returns_empty():
    atom_site = FakeAtomSite(
        columns=["auth_asym_id", "auth_seq_id"],
        rows=[("A", "1")],
    )
    mapping, stats = hook._build_mapping(atom_site)

    assert mapping == {}
    assert "label_asym_id" in stats["missing_fields"]


def test_combined_insertion_and_label_seq():
    atom_site = FakeAtomSite(
        columns=[
            "label_asym_id",
            "auth_asym_id",
            "auth_seq_id",
            "pdbx_PDB_ins_code",
            "label_seq_id",
        ],
        rows=[("A", "H", "100", "A", "52")],
    )
    mapping, _ = hook._build_mapping(atom_site)

    assert mapping[("auth", "H", 100, "A")] == "A"
    assert mapping[("label", "A", 52, "")] == "A"
