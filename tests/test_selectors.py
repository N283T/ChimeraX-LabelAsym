"""Unit tests for dynamic ``la_<label>`` selector registration."""

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


@pytest.fixture
def fake_atomic(monkeypatch):
    """Stub ``chimerax.atomic`` with a synthetic AtomicStructure class."""

    class _AS:
        pass

    fake_mod = MagicMock()
    fake_mod.AtomicStructure = _AS
    monkeypatch.setitem(sys.modules, "chimerax", MagicMock())
    monkeypatch.setitem(sys.modules, "chimerax.atomic", fake_mod)
    return _AS


def test_valid_selector_suffix_accepts_common_labels():
    assert hook._is_valid_selector_suffix("A")
    assert hook._is_valid_selector_suffix("AA")
    assert hook._is_valid_selector_suffix("E1")
    assert hook._is_valid_selector_suffix("chain-1")
    # Digits-only is structurally valid; only the full `la_<label>` name is
    # ever registered, and `la` already satisfies the lead-alpha invariant.
    assert hook._is_valid_selector_suffix("1A")


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

    registered = [
        c.args[0] for c in fake_core_commands.register_selector.call_args_list
    ]
    assert registered == ["la_B"]


def test_register_selectors_skips_invalid_labels(fake_core_commands):
    session = MagicMock()

    hook._register_selectors_for_labels(session, {"A", "bad name", "ok2"})

    registered = sorted(
        c.args[0] for c in fake_core_commands.register_selector.call_args_list
    )
    assert registered == ["la_A", "la_ok2"]


def test_register_selectors_info_logs_skipped_invalid_label_once(fake_core_commands):
    session = MagicMock()

    hook._register_selectors_for_labels(session, {"bad name"})
    hook._register_selectors_for_labels(session, {"bad name"})

    assert session.logger.info.call_count == 1
    assert "bad name" in session.logger.info.call_args.args[0]


def test_register_selectors_warns_on_exception(fake_core_commands):
    session = MagicMock()
    fake_core_commands.register_selector.side_effect = RuntimeError("boom")

    hook._register_selectors_for_labels(session, {"A"})

    session.logger.warning.assert_called_once()
    assert "la_A" in session.logger.warning.call_args.args[0]
    assert "la_A" not in hook._REGISTERED_SELECTORS
    assert "la_A" in hook._FAILED_SELECTORS


def test_register_selectors_does_not_retry_failed(fake_core_commands):
    """A deterministic failure should not re-log on subsequent opens."""
    session = MagicMock()
    fake_core_commands.register_selector.side_effect = RuntimeError("boom")

    hook._register_selectors_for_labels(session, {"A"})
    hook._register_selectors_for_labels(session, {"A"})

    assert fake_core_commands.register_selector.call_count == 1
    assert session.logger.warning.call_count == 1


def test_register_selectors_continues_after_failure(fake_core_commands):
    """One bad label should not prevent sibling labels from registering."""
    session = MagicMock()
    fake_core_commands.register_selector.side_effect = [RuntimeError("boom"), None]

    hook._register_selectors_for_labels(session, ["A", "B"])  # ordered list

    assert fake_core_commands.register_selector.call_count == 2
    assert hook._REGISTERED_SELECTORS == {"la_B"}
    assert hook._FAILED_SELECTORS == {"la_A"}


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


def test_deregister_clears_set_even_when_deregister_raises(fake_core_commands):
    """Uninstall must not leave the set populated even if teardown raises."""
    session = MagicMock()
    fake_core_commands.deregister_selector.side_effect = RuntimeError("boom")
    hook._REGISTERED_SELECTORS.update({"la_A", "la_B"})

    hook._deregister_all_selectors(session)

    assert hook._REGISTERED_SELECTORS == set()
    # Exceptions logged at info level, not silently swallowed.
    assert session.logger.info.call_count == 2


def test_selector_callback_matches_by_label(fake_atomic):
    """The callback returned by ``_make_selector_callback`` filters residues."""
    res_a = MagicMock(label_asym_id="A", atoms="atoms_a")
    res_e = MagicMock(label_asym_id="E", atoms="atoms_e")
    res_unset = MagicMock(label_asym_id=None, atoms="atoms_unset")
    structure = fake_atomic()
    structure.residues = [res_a, res_e, res_unset]

    results = MagicMock()
    callback = hook._make_selector_callback("E")
    callback(MagicMock(), [structure], results)

    results.add_atoms.assert_called_once_with("atoms_e", bonds=True)


def test_selector_callback_skips_non_atomic_structure(fake_atomic):
    """Non-AtomicStructure models (volumes, markers) must be ignored."""

    class _NotStructure:
        residues = "should-not-be-touched"

    non_struct = _NotStructure()
    structure = fake_atomic()
    structure.residues = [MagicMock(label_asym_id="E", atoms="atoms_e")]

    results = MagicMock()
    callback = hook._make_selector_callback("E")
    callback(MagicMock(), [non_struct, structure], results)

    results.add_atoms.assert_called_once_with("atoms_e", bonds=True)


def test_selector_callback_handles_multiple_structures(fake_atomic):
    """Outer loop must visit every structure, not break after the first."""
    s1 = fake_atomic()
    s1.residues = [MagicMock(label_asym_id="E", atoms="atoms_e1")]
    s2 = fake_atomic()
    s2.residues = [MagicMock(label_asym_id="E", atoms="atoms_e2")]

    results = MagicMock()
    callback = hook._make_selector_callback("E")
    callback(MagicMock(), [s1, s2], results)

    assert results.add_atoms.call_count == 2


def test_selector_callback_continues_on_model_failure(fake_atomic):
    """One broken model should not prevent collection from sibling models."""

    class _BrokenStructure(fake_atomic):
        @property
        def residues(self):
            raise RuntimeError("C++ layer gone")

    broken = _BrokenStructure()
    healthy = fake_atomic()
    healthy.residues = [MagicMock(label_asym_id="E", atoms="atoms_good")]

    results = MagicMock()
    session = MagicMock()
    callback = hook._make_selector_callback("E")
    callback(session, [broken, healthy], results)

    results.add_atoms.assert_called_once_with("atoms_good", bonds=True)
    session.logger.warning.assert_called_once()
