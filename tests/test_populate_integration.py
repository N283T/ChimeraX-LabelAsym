"""Integration tests wiring ``_populate`` to selector registration.

Unit tests cover each helper in isolation; this module asserts the end-to-end
flow so a regression where ``_populate`` stops calling
``_register_selectors_for_labels`` (or passes the wrong labels) is caught.
"""

import sys
from unittest.mock import MagicMock

import pytest

import hook


@pytest.fixture(autouse=True)
def reset_state():
    hook._REGISTERED_SELECTORS.clear()
    hook._FAILED_SELECTORS.clear()
    hook._SKIPPED_INVALID_LABELS.clear()
    yield
    hook._REGISTERED_SELECTORS.clear()
    hook._FAILED_SELECTORS.clear()
    hook._SKIPPED_INVALID_LABELS.clear()


class _FakeAtomSite:
    def __init__(self, columns, rows):
        self._columns = list(columns)
        self._rows = rows

    def has_field(self, name):
        return name in self._columns

    def fields(self, requested):
        idx = [self._columns.index(f) for f in requested]
        return [[row[i] for i in idx] for row in self._rows]


def _make_residue(chain_id, number, ins=""):
    r = MagicMock()
    r.chain_id = chain_id
    r.number = number
    r.insertion_code = ins
    # setattr(res, "label_asym_id", ...) mutates the mock normally
    r.label_asym_id = None
    return r


@pytest.fixture
def chimerax_stub(monkeypatch):
    """Stub chimerax.atomic + chimerax.core.commands for _populate."""

    class _AS:  # AtomicStructure stand-in
        pass

    fake_atomic = MagicMock()
    fake_atomic.AtomicStructure = _AS
    fake_commands = MagicMock()
    fake_commands.register_selector = MagicMock()
    fake_commands.deregister_selector = MagicMock()

    monkeypatch.setitem(sys.modules, "chimerax", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.atomic", fake_atomic)
    monkeypatch.setitem(sys.modules, "chimerax.core", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.core.commands", fake_commands)

    return _AS, fake_commands


def _install_fake_extractor(monkeypatch, atom_site):
    """Bypass the mmCIF extraction chain by patching both extractors."""
    monkeypatch.setattr(
        hook, "_extract_atom_site_from_metadata", lambda s, st: atom_site
    )
    monkeypatch.setattr(
        hook, "_extract_atom_site_from_file", lambda s, st: None
    )


def test_populate_registers_selectors_for_observed_labels(
    chimerax_stub, monkeypatch
):
    AS_cls, fake_commands = chimerax_stub
    atom_site = _FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id"],
        rows=[("A", "X", "1"), ("E", "X", "142")],
    )
    _install_fake_extractor(monkeypatch, atom_site)

    structure = AS_cls()
    structure.name = "fake"
    structure.num_residues = 2
    structure.residues = [
        _make_residue("X", 1),
        _make_residue("X", 142),
    ]

    session = MagicMock()
    hook._populate(session, structure)

    registered = sorted(
        c.args[0] for c in fake_commands.register_selector.call_args_list
    )
    assert registered == ["la_A", "la_E"]
    assert hook._REGISTERED_SELECTORS == {"la_A", "la_E"}


def test_populate_skips_selector_registration_when_no_residues_match(
    chimerax_stub, monkeypatch
):
    """mapping has entries, but residues' chain_ids don't match any key."""
    AS_cls, fake_commands = chimerax_stub
    atom_site = _FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id"],
        rows=[("A", "X", "1")],
    )
    _install_fake_extractor(monkeypatch, atom_site)

    structure = AS_cls()
    structure.name = "fake"
    structure.num_residues = 1
    # Residue's chain_id doesn't match the mapping's auth chain "X"
    structure.residues = [_make_residue("OTHER", 1)]

    session = MagicMock()
    hook._populate(session, structure)

    fake_commands.register_selector.assert_not_called()
    assert hook._REGISTERED_SELECTORS == set()


def test_populate_wraps_selector_registration_failure(
    chimerax_stub, monkeypatch
):
    """If register_selector blows up, _populate must still complete the
    attribute assignment and surface a phase-specific warning."""
    AS_cls, fake_commands = chimerax_stub
    atom_site = _FakeAtomSite(
        columns=["label_asym_id", "auth_asym_id", "auth_seq_id"],
        rows=[("A", "X", "1")],
    )
    _install_fake_extractor(monkeypatch, atom_site)

    # Simulate a catastrophic failure that escapes the per-label try/except —
    # e.g. monkey-patch _register_selectors_for_labels itself to raise so the
    # outer phase-specific try/except in _populate is what we're testing.
    def _boom(session, labels):
        raise RuntimeError("mass registration failure")

    monkeypatch.setattr(hook, "_register_selectors_for_labels", _boom)

    residue = _make_residue("X", 1)
    structure = AS_cls()
    structure.name = "fake"
    structure.num_residues = 1
    structure.residues = [residue]

    session = MagicMock()
    hook._populate(session, structure)

    # Attribute still assigned
    assert residue.label_asym_id == "A"
    # Phase-specific warning mentioning selector registration
    warning_msgs = [c.args[0] for c in session.logger.warning.call_args_list]
    assert any("selector registration failed" in m for m in warning_msgs)
