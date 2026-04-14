"""Unit tests for dynamic ``la_<label>`` selector registration."""

import sys
from unittest.mock import MagicMock

import pytest

import hook


@pytest.fixture(autouse=True)
def reset_state():
    hook._REGISTERED_SELECTORS.clear()
    yield
    hook._REGISTERED_SELECTORS.clear()


@pytest.fixture
def fake_core_commands(monkeypatch):
    """Stub ``chimerax.core.commands`` so register/deregister can be captured."""
    fake_mod = MagicMock()
    fake_mod.register_selector = MagicMock()
    fake_mod.deregister_selector = MagicMock()
    monkeypatch.setitem(sys.modules, "chimerax", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.core", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.core.commands", fake_mod)
    return fake_mod


def test_valid_selector_suffix_accepts_common_labels():
    assert hook._is_valid_selector_suffix("A")
    assert hook._is_valid_selector_suffix("AA")
    assert hook._is_valid_selector_suffix("E1")
    assert hook._is_valid_selector_suffix("chain-1")


def test_valid_selector_suffix_rejects_bad_input():
    assert not hook._is_valid_selector_suffix("")
    assert not hook._is_valid_selector_suffix("A B")
    assert not hook._is_valid_selector_suffix("A'")
    assert not hook._is_valid_selector_suffix("A/B")
    assert not hook._is_valid_selector_suffix("A.B")


def test_register_selectors_creates_one_per_unique_label(fake_core_commands):
    session = MagicMock()

    hook._register_selectors_for_labels(session, {"A", "B", "E"})

    names = sorted(
        call.args[0] for call in fake_core_commands.register_selector.call_args_list
    )
    assert names == ["la_A", "la_B", "la_E"]
    assert hook._REGISTERED_SELECTORS == {"la_A", "la_B", "la_E"}


def test_register_selectors_skips_already_registered(fake_core_commands):
    session = MagicMock()
    hook._REGISTERED_SELECTORS.add("la_A")

    hook._register_selectors_for_labels(session, {"A", "B"})

    registered = [c.args[0] for c in fake_core_commands.register_selector.call_args_list]
    assert registered == ["la_B"]


def test_register_selectors_skips_invalid_labels(fake_core_commands):
    session = MagicMock()

    hook._register_selectors_for_labels(session, {"A", "bad name", "ok2"})

    registered = sorted(
        c.args[0] for c in fake_core_commands.register_selector.call_args_list
    )
    assert registered == ["la_A", "la_ok2"]


def test_register_selectors_warns_on_exception(fake_core_commands):
    session = MagicMock()
    fake_core_commands.register_selector.side_effect = RuntimeError("boom")

    hook._register_selectors_for_labels(session, {"A"})

    session.logger.warning.assert_called_once()
    assert "la_A" not in hook._REGISTERED_SELECTORS


def test_deregister_noop_when_empty():
    """No chimerax.core.commands import should be attempted if set is empty."""
    session = MagicMock()
    hook._deregister_all_selectors(session)


def test_deregister_clears_all(fake_core_commands):
    session = MagicMock()
    hook._REGISTERED_SELECTORS.update({"la_A", "la_B"})

    hook._deregister_all_selectors(session)

    assert hook._REGISTERED_SELECTORS == set()
    assert fake_core_commands.deregister_selector.call_count == 2


def test_selector_callback_matches_by_label():
    """The callback returned by ``_make_selector_callback`` filters residues."""
    from unittest.mock import MagicMock as Mock

    fake_atomic = MagicMock()

    class _AS:  # stands in for chimerax.atomic.AtomicStructure
        pass

    fake_atomic.AtomicStructure = _AS
    sys.modules.setdefault("chimerax", MagicMock())
    sys.modules["chimerax.atomic"] = fake_atomic

    res_a = Mock(label_asym_id="A", atoms="atoms_a")
    res_e = Mock(label_asym_id="E", atoms="atoms_e")
    res_unset = Mock(label_asym_id=None, atoms="atoms_unset")
    structure = _AS()
    structure.residues = [res_a, res_e, res_unset]

    results = Mock()
    callback = hook._make_selector_callback("E")
    callback(MagicMock(), [structure], results)

    results.add_atoms.assert_called_once_with("atoms_e", bonds=True)
