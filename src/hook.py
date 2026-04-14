"""Populate a ``label_asym_id`` attribute on residues from mmCIF metadata.

Strategy:
1. On install, register ``Residue.label_asym_id`` as a persistent string
   attribute and subscribe to the ``ADD_MODELS`` trigger.
2. For each new ``AtomicStructure``, pull the ``atom_site`` CIF table from
   the structure metadata (``get_mmcif_tables_from_metadata``) or, if that
   fails, re-read the source text mmCIF via ``get_cif_tables``.
3. Build a lookup keyed by ``(scheme, chain, seq_num, ins_code) -> label``
   for both the auth scheme and the label scheme. Residues are then matched
   regardless of whether the file was opened with ``prefer_auth`` true or
   false (ChimeraX's ``Residue.chain_id`` / ``number`` follow that flag).
4. Assign ``label_asym_id`` to every residue that matches either key.
   Diagnostics are logged: the number of matched residues, skipped rows
   with non-integer sequence ids, and any structure that carried
   ``atom_site`` data but produced zero residue matches (likely a key-
   scheme mismatch). Non-mmCIF structures are silently skipped.
5. Register a ``la_<label>`` selector via
   ``chimerax.core.commands.register_selector`` for each unique label
   observed on at least one residue (filtered by ``_is_valid_selector_suffix``
   so unusual mmCIF ids cannot produce invalid selector names). Selectors
   are session-global: opening multiple structures that share a label id
   reuses the existing selector, and ``select la_X`` matches residues
   across every open model. ``uninstall()`` deregisters every selector
   this module registered.
"""

from __future__ import annotations

from pathlib import Path

ATTR_NAME = "label_asym_id"
SELECTOR_PREFIX = "la_"
_TEXT_MMCIF_SUFFIXES = frozenset({".cif", ".mmcif"})
_REGISTERED = False
_HANDLER = None
_MMCIF_WARNING_LOGGED = False
# Selectors successfully registered with ChimeraX. Session-global by design:
# multiple structures commonly share labels, so deregistering on model close
# would break selectors for still-open structures. Cleared only on uninstall.
_REGISTERED_SELECTORS: set[str] = set()
# Labels whose registration raised — tracked to avoid spamming the log on
# every subsequent open when the failure is deterministic (e.g. name clash).
_FAILED_SELECTORS: set[str] = set()
# Labels rejected by _is_valid_selector_suffix; deduped so we only info-log
# once per skipped label per session.
_SKIPPED_INVALID_LABELS: set[str] = set()


def _log_mmcif_missing_once(session) -> None:
    global _MMCIF_WARNING_LOGGED
    if _MMCIF_WARNING_LOGGED:
        return
    session.logger.warning(
        "[label-asym] chimerax.mmcif is not importable; "
        "label_asym_id attributes will not be populated"
    )
    _MMCIF_WARNING_LOGGED = True


def _register_attr_once(session) -> None:
    global _REGISTERED
    if _REGISTERED:
        return
    from chimerax.atomic import Residue

    Residue.register_attr(session, ATTR_NAME, "ChimeraX-LabelAsym", attr_type=str)
    _REGISTERED = True


def _extract_atom_site_from_metadata(session, structure):
    """Return the atom_site CIFTable from structure metadata, or None."""
    try:
        from chimerax.mmcif import get_mmcif_tables_from_metadata
    except ImportError:
        _log_mmcif_missing_once(session)
        return None
    try:
        tables = get_mmcif_tables_from_metadata(structure, ["atom_site"])
    except (KeyError, ValueError, AttributeError) as exc:
        session.logger.info(
            f"[label-asym] {structure.name}: no atom_site in metadata ({exc})"
        )
        return None
    if not tables:
        return None
    return tables[0]


def _extract_atom_site_from_file(session, structure):
    """Fallback: re-parse the source text mmCIF from disk."""
    filename = getattr(structure, "filename", None)
    if not filename:
        return None
    path = Path(filename)
    if not path.is_file():
        return None
    if path.suffix.lower() not in _TEXT_MMCIF_SUFFIXES:
        return None
    try:
        from chimerax.mmcif import get_cif_tables
    except ImportError:
        _log_mmcif_missing_once(session)
        return None
    try:
        tables = get_cif_tables(str(path), ["atom_site"])
    except OSError as exc:
        session.logger.info(f"[label-asym] cannot re-read {path.name}: {exc}")
        return None
    except (KeyError, ValueError) as exc:
        session.logger.warning(f"[label-asym] parse failed for {path.name}: {exc}")
        return None
    if not tables:
        return None
    return tables[0]


def _build_mapping(atom_site) -> tuple[dict, dict]:
    """Build a key lookup from an ``atom_site`` CIFTable.

    Returns ``(mapping, stats)`` where ``mapping`` has entries keyed by
    ``(scheme, chain, seq_num, ins_code) -> label_asym_id`` and ``stats``
    reports the number of rows skipped due to non-numeric sequence ids.
    Scheme is ``"auth"`` for the auth_asym_id + auth_seq_id key and
    ``"label"`` for the label_asym_id + label_seq_id key. Both are emitted
    so we can match residues regardless of ``prefer_auth``.
    """
    stats = {"skipped_auth_seq": 0, "skipped_label_seq": 0, "missing_fields": []}

    required = ("label_asym_id", "auth_asym_id", "auth_seq_id")
    missing = [f for f in required if not atom_site.has_field(f)]
    if missing:
        stats["missing_fields"] = missing
        return {}, stats

    want = list(required)
    has_ins = atom_site.has_field("pdbx_PDB_ins_code")
    has_label_seq = atom_site.has_field("label_seq_id")
    if has_ins:
        want.append("pdbx_PDB_ins_code")
    if has_label_seq:
        want.append("label_seq_id")

    mapping: dict = {}
    for row in atom_site.fields(want):
        idx = 0
        lbl = row[idx]
        idx += 1
        auth = row[idx]
        idx += 1
        aseq = row[idx]
        idx += 1
        ins = ""
        if has_ins:
            ins_raw = row[idx]
            idx += 1
            ins = ins_raw if ins_raw and ins_raw not in ("?", ".") else ""

        try:
            aseq_num = int(aseq)
            mapping.setdefault(("auth", auth, aseq_num, ins), lbl)
        except (TypeError, ValueError):
            stats["skipped_auth_seq"] += 1

        if has_label_seq:
            lseq = row[idx]
            idx += 1
            try:
                lseq_num = int(lseq)
                mapping.setdefault(("label", lbl, lseq_num, ""), lbl)
            except (TypeError, ValueError):
                stats["skipped_label_seq"] += 1

    return mapping, stats


def _is_valid_selector_suffix(label: str) -> bool:
    """Check whether ``la_<label>`` would pass ChimeraX's selector name check.

    ChimeraX (``chimerax.core.commands.atomspec.register_selector``) requires
    ``name[0].isalpha()`` and every subsequent character to satisfy
    ``c.isalnum() or c in "-+_"``. The fixed ``la_`` prefix satisfies the lead
    and tail invariants; this function validates only the user-supplied
    ``<label>`` portion, so unusual mmCIF ids (``HOH 1``, ``A'``) are caught
    before we emit a ``register_selector`` warning.
    """
    if not label:
        return False
    return all(c.isalnum() or c in "-+_" for c in label)


def _make_selector_callback(label: str):
    """Build a selector callback for a single ``label_asym_id`` value.

    The returned closure has the signature required by
    ``chimerax.core.commands.register_selector``: it receives the session, a
    list of candidate ``Model`` instances, and a mutable ``Objects``
    accumulator. It walks each ``AtomicStructure`` model's residues and adds
    the atoms (plus bonds, so ribbons/sticks follow) of any residue whose
    ``label_asym_id`` attribute equals ``label``.

    Failures in one model (e.g. partially-torn-down structure with a bad
    ``residues`` accessor) are logged and do not abort the remaining models,
    so partial results still land in ``results``.
    """

    def _select(session, models, results):
        import traceback

        from chimerax.atomic import AtomicStructure

        for m in models:
            if not isinstance(m, AtomicStructure):
                continue
            try:
                for res in m.residues:
                    if getattr(res, ATTR_NAME, None) == label:
                        results.add_atoms(res.atoms, bonds=True)
            except Exception:
                tb = traceback.format_exc()
                session.logger.warning(
                    f"[label-asym] selector la_{label} failed on "
                    f"{getattr(m, 'name', m)}; continuing:\n{tb}"
                )

    return _select


def _register_selectors_for_labels(session, labels) -> None:
    """Register a ``la_<label>`` selector for each label in ``labels``.

    Idempotent across calls: ``_REGISTERED_SELECTORS`` deduplicates so that
    re-opening the same structure (or opening a second structure that shares
    labels) is a no-op. Labels whose characters would be rejected by
    ChimeraX's ``register_selector`` lead/tail checks are skipped via
    ``_is_valid_selector_suffix`` and info-logged once via
    ``_SKIPPED_INVALID_LABELS``. Labels whose registration raised are tracked
    in ``_FAILED_SELECTORS`` so we do not re-spam the log on subsequent opens
    when the failure is deterministic (e.g. selector-name collision).
    """
    import traceback

    from chimerax.core.commands import register_selector

    for lbl in labels:
        if not _is_valid_selector_suffix(lbl):
            if lbl not in _SKIPPED_INVALID_LABELS:
                _SKIPPED_INVALID_LABELS.add(lbl)
                session.logger.info(
                    f"[label-asym] label {lbl!r} is not a valid selector "
                    f"suffix (contains whitespace or punctuation); "
                    f'use `select ::label_asym_id="{lbl}"` instead'
                )
            continue
        name = f"{SELECTOR_PREFIX}{lbl}"
        if name in _REGISTERED_SELECTORS or name in _FAILED_SELECTORS:
            continue
        try:
            register_selector(
                name,
                _make_selector_callback(lbl),
                session.logger,
                desc=f"residues with label_asym_id={lbl}",
            )
        except Exception:
            tb = traceback.format_exc()
            session.logger.warning(
                f"[label-asym] failed to register selector {name}; "
                f'fall back to `select ::label_asym_id="{lbl}"`:\n{tb}'
            )
            _FAILED_SELECTORS.add(name)
            continue
        _REGISTERED_SELECTORS.add(name)


def _deregister_all_selectors(session) -> None:
    """Deregister every selector this module registered.

    Runs during bundle teardown (``uninstall``). Exceptions are logged at
    ``info`` level rather than swallowed silently — shutdown paths are
    best-effort but not voiceless. ``_REGISTERED_SELECTORS`` is cleared
    unconditionally so a subsequent ``install`` in the same process starts
    from a clean slate; ``_FAILED_SELECTORS`` and ``_SKIPPED_INVALID_LABELS``
    are reset too so previously-failed names can be retried.
    """
    if not _REGISTERED_SELECTORS:
        _FAILED_SELECTORS.clear()
        _SKIPPED_INVALID_LABELS.clear()
        return
    import traceback

    from chimerax.core.commands import deregister_selector

    for name in list(_REGISTERED_SELECTORS):
        try:
            deregister_selector(name, session.logger)
        except Exception:
            tb = traceback.format_exc()
            session.logger.info(
                f"[label-asym] deregister_selector({name}) failed during "
                f"shutdown:\n{tb}"
            )
    _REGISTERED_SELECTORS.clear()
    _FAILED_SELECTORS.clear()
    _SKIPPED_INVALID_LABELS.clear()


def _populate(session, structure) -> None:
    from chimerax.atomic import AtomicStructure

    if not isinstance(structure, AtomicStructure):
        return

    atom_site = _extract_atom_site_from_metadata(session, structure)
    if atom_site is None:
        atom_site = _extract_atom_site_from_file(session, structure)
    if atom_site is None:
        return

    mapping, stats = _build_mapping(atom_site)

    if stats["missing_fields"]:
        session.logger.info(
            f"[label-asym] {structure.name}: atom_site lacks "
            f"{', '.join(stats['missing_fields'])}; skipping"
        )
        return

    if not mapping:
        session.logger.info(
            f"[label-asym] {structure.name}: atom_site yielded no usable rows"
        )
        return

    hit = 0
    seen_labels: set[str] = set()
    for res in structure.residues:
        ins = res.insertion_code or ""
        lbl = mapping.get(("auth", res.chain_id, res.number, ins))
        if lbl is None:
            lbl = mapping.get(("label", res.chain_id, res.number, ins))
        if lbl is not None:
            setattr(res, ATTR_NAME, lbl)
            seen_labels.add(lbl)
            hit += 1

    if seen_labels:
        try:
            _register_selectors_for_labels(session, seen_labels)
        except Exception:
            import traceback

            tb = traceback.format_exc()
            session.logger.warning(
                f"[label-asym] {structure.name}: attribute assigned but "
                f"selector registration failed:\n{tb}"
            )

    total = structure.num_residues
    if stats["skipped_auth_seq"]:
        session.logger.info(
            f"[label-asym] {structure.name}: skipped "
            f"{stats['skipped_auth_seq']} atom_site rows with non-numeric "
            "auth_seq_id"
        )

    if hit == 0:
        session.logger.warning(
            f"[label-asym] {structure.name}: atom_site present but 0 of "
            f"{total} residues matched; check prefer_auth and ins_code handling"
        )
    else:
        session.logger.info(
            f"[label-asym] {structure.name}: assigned label_asym_id to "
            f"{hit}/{total} residues"
        )


def _on_add_models(session, trigger_name, models):
    import traceback

    for m in models:
        try:
            _populate(session, m)
        except Exception:
            tb = traceback.format_exc()
            session.logger.warning(
                f"[label-asym] failed to populate {getattr(m, 'name', m)}:\n{tb}"
            )


def install(session) -> None:
    _register_attr_once(session)
    global _HANDLER
    if _HANDLER is not None:
        return
    from chimerax.core.models import ADD_MODELS

    _HANDLER = session.triggers.add_handler(
        ADD_MODELS, lambda tname, mlist: _on_add_models(session, tname, mlist)
    )


def uninstall(session) -> None:
    global _HANDLER
    _deregister_all_selectors(session)
    if _HANDLER is None:
        return
    try:
        session.triggers.remove_handler(_HANDLER)
    except (ValueError, KeyError) as exc:
        session.logger.warning(f"[label-asym] handler removal failed: {exc}")
    finally:
        _HANDLER = None
