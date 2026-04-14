"""Unit tests for ``install`` / ``uninstall`` lifecycle."""

import sys
from unittest.mock import MagicMock

import pytest

import hook


@pytest.fixture(autouse=True)
def reset_state():
    """Reset module-level state between tests so they do not leak."""
    hook._REGISTERED = False
    hook._HANDLER = None
    hook._MMCIF_WARNING_LOGGED = False
    yield
    hook._HANDLER = None


def _fake_chimerax_modules(monkeypatch):
    """Stub out the ChimeraX imports that ``install`` pulls in."""
    fake_chimerax = MagicMock()
    fake_core = MagicMock()
    fake_atomic = MagicMock()
    fake_models = MagicMock()
    fake_models.ADD_MODELS = "add models"
    monkeypatch.setitem(sys.modules, "chimerax", fake_chimerax)
    monkeypatch.setitem(sys.modules, "chimerax.core", fake_core)
    monkeypatch.setitem(sys.modules, "chimerax.core.models", fake_models)
    monkeypatch.setitem(sys.modules, "chimerax.atomic", fake_atomic)
    return fake_atomic, fake_models


def _fake_session():
    session = MagicMock()
    session.triggers.add_handler.return_value = "handler-1"
    return session


def test_install_subscribes_once(monkeypatch):
    _fake_chimerax_modules(monkeypatch)
    session = _fake_session()

    hook.install(session)
    hook.install(session)

    assert session.triggers.add_handler.call_count == 1
    assert hook._HANDLER == "handler-1"


def test_uninstall_removes_registered_handler(monkeypatch):
    _fake_chimerax_modules(monkeypatch)
    session = _fake_session()

    hook.install(session)
    hook.uninstall(session)

    session.triggers.remove_handler.assert_called_once_with("handler-1")
    assert hook._HANDLER is None


def test_uninstall_without_install_is_noop():
    session = MagicMock()
    hook.uninstall(session)

    session.triggers.remove_handler.assert_not_called()


def test_uninstall_logs_when_remove_handler_raises(monkeypatch):
    _fake_chimerax_modules(monkeypatch)
    session = _fake_session()
    hook.install(session)

    session.triggers.remove_handler.side_effect = ValueError("already gone")
    hook.uninstall(session)

    session.logger.warning.assert_called_once()
    assert hook._HANDLER is None
