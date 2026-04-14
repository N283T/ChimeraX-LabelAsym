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


def _make_atoms(labels):
    """Build a minimal Atoms-like mock with a residues list and colors array."""
    residues = [MagicMock(label_asym_id=lbl) for lbl in labels]
    atoms = MagicMock()
    atoms.__len__ = lambda self: len(labels)
    atoms.residues = residues
    atoms.colors = np.full((len(labels), 4), 255, dtype=np.uint8)

    def _filter(mask):
        idxs = np.where(mask)[0]
        sub = MagicMock()
        sub.__len__ = lambda self: len(idxs)
        sub.colors = atoms.colors[idxs].copy()

        def _set_colors(new):
            atoms.colors[idxs] = new

        type(sub).colors = property(
            lambda self: atoms.colors[idxs].copy(),
            lambda self, v: _set_colors(v),
        )
        return sub

    atoms.filter = _filter
    return atoms


def test_labelcolor_assigns_distinct_rgb_per_label(stubbed_chimerax):
    import commands as commands_mod

    session = MagicMock()
    atoms = _make_atoms(["A", "A", "E", "E"])
    objects = MagicMock(atoms=atoms)

    commands_mod.labelcolor(session, objects)

    assert not (atoms.colors[0, :3] == atoms.colors[2, :3]).all()
    assert (atoms.colors[0, :3] == atoms.colors[1, :3]).all()
    assert (atoms.colors[2, :3] == atoms.colors[3, :3]).all()


def test_labelcolor_skips_atoms_without_label(stubbed_chimerax):
    import commands as commands_mod

    session = MagicMock()
    atoms = _make_atoms(["A", "", None, "E"])
    objects = MagicMock(atoms=atoms)

    commands_mod.labelcolor(session, objects)

    # unlabeled atoms keep default (255) in RGB; labeled atoms are recolored
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
    atoms = _make_atoms(["", "", None])
    objects = MagicMock(atoms=atoms)

    commands_mod.labelcolor(session, objects)

    session.logger.warning.assert_called_once()


def test_labelcolor_defaults_objects_to_all_objects(stubbed_chimerax):
    _, fake_objects = stubbed_chimerax
    import commands as commands_mod

    atoms = _make_atoms(["A"])
    fake_objects.all_objects.return_value = MagicMock(atoms=atoms)

    session = MagicMock()
    commands_mod.labelcolor(session)

    fake_objects.all_objects.assert_called_once_with(session)
