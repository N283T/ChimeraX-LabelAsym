"""Unit tests for ``_extract_atom_site_from_file`` suffix gating."""

from unittest.mock import MagicMock

import hook


class _FakeStructure:
    def __init__(self, filename=None, name="fake"):
        self.filename = filename
        self.name = name


def test_no_filename_returns_none():
    session = MagicMock()
    assert hook._extract_atom_site_from_file(session, _FakeStructure(None)) is None


def test_missing_file_returns_none():
    session = MagicMock()
    assert (
        hook._extract_atom_site_from_file(session, _FakeStructure("/nonexistent.cif"))
        is None
    )


def test_pdb_suffix_is_skipped(tmp_path):
    p = tmp_path.joinpath("test.pdb")
    p.write_text("HEADER\n")
    session = MagicMock()

    assert hook._extract_atom_site_from_file(session, _FakeStructure(str(p))) is None


def test_bcif_suffix_is_skipped(tmp_path):
    """bCIF must not be routed through the text-only ``get_cif_tables`` path."""
    p = tmp_path.joinpath("test.bcif")
    p.write_bytes(b"\x00\x01")
    session = MagicMock()

    assert hook._extract_atom_site_from_file(session, _FakeStructure(str(p))) is None


def test_cif_suffix_is_accepted_for_parsing(tmp_path, monkeypatch):
    """A .cif path proceeds past suffix gating and calls ``get_cif_tables``."""
    import sys

    p = tmp_path.joinpath("test.cif")
    p.write_text("data_test\n")
    session = MagicMock()

    sentinel = object()
    fake_mmcif = MagicMock(get_cif_tables=MagicMock(return_value=[sentinel]))
    monkeypatch.setitem(sys.modules, "chimerax", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.mmcif", fake_mmcif)

    result = hook._extract_atom_site_from_file(session, _FakeStructure(str(p)))

    assert result is sentinel
    fake_mmcif.get_cif_tables.assert_called_once_with(str(p), ["atom_site"])
