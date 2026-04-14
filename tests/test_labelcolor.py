"""Unit tests for the ``labelcolor`` command."""

import sys
from unittest.mock import MagicMock

import numpy as np
import pytest


@pytest.fixture
def stubbed_chimerax(monkeypatch):
    """Stub chimerax modules imported by ``labelcolor`` at call time."""
    fake_colors = MagicMock()

    def _chain_colors(ids):
        """Deterministic palette: one unique RGBA per unique id."""
        seen: dict[str, tuple[int, int, int, int]] = {}
        rgba = np.zeros((len(ids), 4), dtype=np.uint8)
        for i, cid in enumerate(ids):
            if cid not in seen:
                n = len(seen) + 1
                seen[cid] = (n * 10, n * 20, n * 30, 255)
            rgba[i] = seen[cid]
        return rgba

    fake_colors.chain_colors = _chain_colors

    fake_objects = MagicMock()
    fake_objects.all_objects = MagicMock(return_value=MagicMock())

    monkeypatch.setitem(sys.modules, "chimerax", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.atomic", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.atomic.colors", fake_colors)
    monkeypatch.setitem(sys.modules, "chimerax.core", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.core.commands", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.core.objects", fake_objects)

    return fake_colors, fake_objects


class _ColorCarrier:
    """Minimal atoms/residues stand-in with readable/writable `.colors` arrays.

    Mirrors ChimeraX's scatter semantics: the setter captures the assigned
    array so tests can assert against it. Uses a plain Python attribute so
    `__len__` works and `MagicMock`'s auto-spec doesn't get in the way.
    """

    def __init__(self, labels, initial_alpha=255):
        self.residues = [MagicMock(label_asym_id=lbl) for lbl in labels]
        self._labels = list(labels)
        rgba = np.full((len(labels), 4), 255, dtype=np.uint8)
        rgba[:, 3] = initial_alpha
        self.colors = rgba

    def __len__(self):
        return len(self._labels)

    def __iter__(self):
        # `_labels_of(residues)` iterates the residues directly when called
        # with a residues collection. Atoms path uses `.residues` first.
        return iter(self.residues)


def _residues_carrier(labels, initial_alpha=200):
    carrier = _ColorCarrier(labels, initial_alpha=initial_alpha)
    rgba = np.full((len(labels), 4), 128, dtype=np.uint8)
    rgba[:, 3] = initial_alpha
    carrier.ribbon_colors = rgba.copy()
    carrier.ring_colors = rgba.copy()
    return carrier


def test_labelcolor_assigns_distinct_rgb_per_label(stubbed_chimerax):
    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["A", "A", "E", "E"])
    objects = MagicMock(atoms=atoms, residues=_residues_carrier([]))

    commands_mod.labelcolor(session, objects)

    assert not (atoms.colors[0, :3] == atoms.colors[2, :3]).all()
    assert (atoms.colors[0, :3] == atoms.colors[1, :3]).all()
    assert (atoms.colors[2, :3] == atoms.colors[3, :3]).all()


def test_labelcolor_preserves_alpha_on_atoms(stubbed_chimerax):
    """Recoloring must leave the alpha column untouched for every atom."""
    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["A", "E"], initial_alpha=128)
    atoms.colors[:, 3] = np.array([128, 200], dtype=np.uint8)
    objects = MagicMock(atoms=atoms, residues=_residues_carrier([]))

    commands_mod.labelcolor(session, objects)

    assert atoms.colors[0, 3] == 128
    assert atoms.colors[1, 3] == 200


def test_labelcolor_writes_ribbon_and_ring_colors(stubbed_chimerax):
    """The residue-side ribbon/ring colors must be recolored too."""
    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["A", "E"])
    residues = _residues_carrier(["A", "E"])
    objects = MagicMock(atoms=atoms, residues=residues)

    commands_mod.labelcolor(session, objects)

    # Ribbon colors got a palette entry per residue (different per label)
    assert not (residues.ribbon_colors[0, :3] == residues.ribbon_colors[1, :3]).all()
    assert not (residues.ring_colors[0, :3] == residues.ring_colors[1, :3]).all()
    # Alpha preserved on ribbons/rings
    assert (residues.ribbon_colors[:, 3] == 200).all()
    assert (residues.ring_colors[:, 3] == 200).all()


def test_labelcolor_skips_ribbon_branch_when_residues_empty(stubbed_chimerax):
    """No crash when the selection has no residues of its own."""
    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["A"])
    # Empty residues collection — len() == 0, branch must be skipped
    empty = _residues_carrier([])
    objects = MagicMock(atoms=atoms, residues=empty)

    commands_mod.labelcolor(session, objects)

    # Atoms still colored
    assert not (atoms.colors[0, :3] == 255).all()


def test_labelcolor_skips_atoms_without_label(stubbed_chimerax):
    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["A", "", None, "E"])
    objects = MagicMock(atoms=atoms, residues=_residues_carrier([]))

    commands_mod.labelcolor(session, objects)

    # Unlabeled atoms keep default RGB (255); labeled atoms are recolored
    assert (atoms.colors[1, :3] == 255).all()
    assert (atoms.colors[2, :3] == 255).all()
    assert not (atoms.colors[0, :3] == 255).all()
    assert not (atoms.colors[3, :3] == 255).all()


def test_labelcolor_warns_when_no_atoms(stubbed_chimerax):
    import commands as commands_mod

    session = MagicMock()
    atoms = MagicMock()
    atoms.__len__ = lambda self: 0
    objects = MagicMock(atoms=atoms)

    commands_mod.labelcolor(session, objects)

    session.logger.warning.assert_called_once()


def test_labelcolor_warns_when_no_labels(stubbed_chimerax):
    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["", "", None])
    objects = MagicMock(atoms=atoms, residues=_residues_carrier([]))

    commands_mod.labelcolor(session, objects)

    session.logger.warning.assert_called_once()


def test_labelcolor_defaults_objects_to_all_objects(stubbed_chimerax):
    _, fake_objects = stubbed_chimerax
    import commands as commands_mod

    atoms = _ColorCarrier(["A"])
    fake_objects.all_objects.return_value = MagicMock(
        atoms=atoms, residues=_residues_carrier([])
    )

    session = MagicMock()
    commands_mod.labelcolor(session)

    fake_objects.all_objects.assert_called_once_with(session)


def test_labelcolor_handles_palette_shape_mismatch(stubbed_chimerax):
    """A broken palette should warn and leave colors unchanged."""
    fake_colors, _ = stubbed_chimerax
    fake_colors.chain_colors = lambda ids: np.zeros((0, 4), dtype=np.uint8)

    import commands as commands_mod

    session = MagicMock()
    atoms = _ColorCarrier(["A", "E"])
    original = atoms.colors.copy()
    objects = MagicMock(atoms=atoms, residues=_residues_carrier([]))

    commands_mod.labelcolor(session, objects)

    session.logger.warning.assert_called_once()
    assert (atoms.colors == original).all()


def test_apply_rgb_rejects_palette_size_mismatch():
    """``_apply_rgb`` should raise a clear ValueError, not numpy broadcast."""
    import commands as commands_mod

    labels = ["A", "B"]
    mask = np.array([True, True])
    current = np.full((2, 4), 255, dtype=np.uint8)

    with pytest.raises(ValueError, match="1 colors for 2 labels"):
        commands_mod._apply_rgb(
            labels,
            mask,
            lambda ids: np.zeros((1, 4), dtype=np.uint8),
            current,
        )
