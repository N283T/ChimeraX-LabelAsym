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
"""

from __future__ import annotations

from pathlib import Path

ATTR_NAME = "label_asym_id"
_TEXT_MMCIF_SUFFIXES = frozenset({".cif", ".mmcif"})
_REGISTERED = False
_HANDLER = None
_MMCIF_WARNING_LOGGED = False


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
    for res in structure.residues:
        ins = res.insertion_code or ""
        lbl = mapping.get(("auth", res.chain_id, res.number, ins))
        if lbl is None:
            lbl = mapping.get(("label", res.chain_id, res.number, ins))
        if lbl is not None:
            setattr(res, ATTR_NAME, lbl)
            hit += 1

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
    if _HANDLER is None:
        return
    try:
        session.triggers.remove_handler(_HANDLER)
    except (ValueError, KeyError) as exc:
        session.logger.warning(f"[label-asym] handler removal failed: {exc}")
    finally:
        _HANDLER = None
