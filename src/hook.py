"""Populate a ``label_asym_id`` attribute on residues from mmCIF metadata.

Strategy:
1. Register ``Residue.label_asym_id`` as a persistent string attribute once.
2. On every ADD_MODELS trigger, walk new AtomicStructures and try to extract
   ``atom_site`` from their mmCIF metadata (``get_mmcif_tables_from_metadata``).
3. If metadata lacks atom_site, fall back to re-reading the source file via
   ``get_cif_tables`` when a plausible mmCIF path is available.
4. Build a ``(auth_asym_id, auth_seq_id, ins_code) -> label_asym_id`` map and
   assign it to each residue. Silent on non-mmCIF structures.
"""

from __future__ import annotations

from pathlib import Path

ATTR_NAME = "label_asym_id"
_REGISTERED = False
_HANDLER = None


def _register_attr_once(session) -> None:
    global _REGISTERED
    if _REGISTERED:
        return
    from chimerax.atomic import Residue

    Residue.register_attr(
        session, ATTR_NAME, "ChimeraX-LabelAsym", attr_type=str
    )
    _REGISTERED = True


def _extract_atom_site_from_metadata(structure):
    """Return the atom_site CIFTable or None."""
    try:
        from chimerax.mmcif import get_mmcif_tables_from_metadata
    except ImportError:
        return None
    try:
        tables = get_mmcif_tables_from_metadata(structure, ["atom_site"])
    except Exception:
        return None
    if not tables:
        return None
    return tables[0]


def _extract_atom_site_from_file(structure):
    """Fallback: re-parse the source mmCIF from disk."""
    filename = getattr(structure, "filename", None)
    if not filename:
        return None
    path = Path(filename)
    if not path.is_file():
        return None
    suffix = path.suffix.lower()
    if suffix not in (".cif", ".mmcif", ".bcif"):
        return None
    try:
        from chimerax.mmcif import get_cif_tables
    except ImportError:
        return None
    try:
        tables = get_cif_tables(str(path), ["atom_site"])
    except Exception:
        return None
    if not tables:
        return None
    return tables[0]


def _build_mapping(atom_site) -> dict:
    """Build (auth_asym_id, auth_seq_id, ins_code) -> label_asym_id."""
    if not atom_site.has_field("label_asym_id"):
        return {}
    if not atom_site.has_field("auth_asym_id") or not atom_site.has_field("auth_seq_id"):
        return {}

    fields = ["label_asym_id", "auth_asym_id", "auth_seq_id"]
    has_ins = atom_site.has_field("pdbx_PDB_ins_code")
    if has_ins:
        fields.append("pdbx_PDB_ins_code")

    mapping: dict = {}
    for row in atom_site.fields(fields):
        if has_ins:
            lbl, auth, seq, ins = row
        else:
            lbl, auth, seq = row
            ins = ""
        try:
            seq_num = int(seq)
        except (TypeError, ValueError):
            continue
        ins_code = ins if ins and ins not in ("?", ".") else ""
        key = (auth, seq_num, ins_code)
        mapping.setdefault(key, lbl)
    return mapping


def _populate(session, structure) -> None:
    from chimerax.atomic import AtomicStructure

    if not isinstance(structure, AtomicStructure):
        return

    atom_site = _extract_atom_site_from_metadata(structure)
    if atom_site is None:
        atom_site = _extract_atom_site_from_file(structure)
    if atom_site is None:
        return

    mapping = _build_mapping(atom_site)
    if not mapping:
        return

    hit = 0
    for res in structure.residues:
        key = (res.chain_id, res.number, res.insertion_code or "")
        lbl = mapping.get(key)
        if lbl is not None:
            setattr(res, ATTR_NAME, lbl)
            hit += 1

    if hit:
        session.logger.info(
            f"[label-asym] {structure.name}: assigned label_asym_id to "
            f"{hit}/{structure.num_residues} residues"
        )


def _on_add_models(session, trigger_name, models):
    for m in models:
        try:
            _populate(session, m)
        except Exception as exc:
            session.logger.warning(
                f"[label-asym] skipped {getattr(m, 'name', m)}: {exc}"
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
    except Exception:
        pass
    _HANDLER = None
