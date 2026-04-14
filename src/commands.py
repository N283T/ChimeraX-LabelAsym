"""ChimeraX commands for the LabelAsym bundle.

Currently provides ``colorbylabel``, a ``color bychain`` analogue that
groups residues by ``label_asym_id`` instead of ``auth_asym_id``. The
default palette matches ``chimerax.atomic.colors.chain_colors`` so heme
chains (label E/F/G/H in 4hhb) get visibly distinct colors from their
associated proteins (label A/B/C/D) while preserving the same palette
semantics users already know from ``color bychain``.

Named ``colorbylabel`` (not ``labelcolor``) because ChimeraX already has
a first-class concept of 3D text "labels" — the name avoids a collision
with that mental model.
"""

from __future__ import annotations

try:
    from .hook import ATTR_NAME
except ImportError as pkg_err:  # tests import bare modules via sys.path
    try:
        from hook import ATTR_NAME
    except ImportError as bare_err:
        raise ImportError(
            f"label-asym: could not import hook as package "
            f"({pkg_err}) or as top-level module ({bare_err})"
        ) from bare_err


def _labels_of(items) -> list[str]:
    """Read ``label_asym_id`` per item, using "" when the attribute is unset."""
    return [getattr(it, ATTR_NAME, None) or "" for it in items]


def _apply_rgb(labels, mask, palette_fn, current_rgba):
    """Return a copy of ``current_rgba`` recolored for ``mask==True`` rows.

    RGB columns (0..2) of the ``mask==True`` rows are replaced by the output
    of ``palette_fn(keep_labels)``. Alpha (column 3) is preserved for every
    row; ``mask==False`` rows are copied unchanged. ``palette_fn`` must
    return an ``(n, 4)`` uint8 array where ``n == sum(mask)``; a shape
    mismatch raises ``ValueError`` so callers can log a specific diagnostic
    instead of a bare numpy broadcast error.
    """
    import numpy as np

    keep = [lbl for lbl, m in zip(labels, mask) if m]
    palette = palette_fn(keep)
    if palette.shape[0] != len(keep):
        raise ValueError(
            f"palette returned {palette.shape[0]} colors for {len(keep)} labels"
        )
    out = current_rgba.copy()
    idx = np.where(mask)[0]
    out[idx, :3] = palette[:, :3]
    return out


def colorbylabel(session, objects=None):
    """Color atoms, ribbons, and rings by ``label_asym_id``."""
    import numpy as np

    from chimerax.atomic.colors import chain_colors
    from chimerax.core.objects import all_objects

    if objects is None:
        objects = all_objects(session)

    atoms = objects.atoms
    if len(atoms) == 0:
        session.logger.warning(
            "[label-asym] colorbylabel: selection contains no atoms; "
            "check your atom specifier or open a structure first"
        )
        return

    atom_labels = _labels_of(atoms.residues)
    atom_mask = np.array([bool(lbl) for lbl in atom_labels], dtype=bool)
    if not atom_mask.any():
        session.logger.warning(
            f"[label-asym] colorbylabel: none of {len(atoms)} atoms have a "
            f"label_asym_id attribute. Either the structure was not opened "
            f"from mmCIF, or atom_site metadata was missing — see earlier "
            f"[label-asym] log messages for details."
        )
        return

    try:
        atoms.colors = _apply_rgb(atom_labels, atom_mask, chain_colors, atoms.colors)
    except (ValueError, IndexError) as exc:
        session.logger.warning(
            f"[label-asym] colorbylabel: palette application failed ({exc}); "
            f"atom colors unchanged"
        )
        return

    residues = objects.residues
    if len(residues) > 0:
        res_labels = _labels_of(residues)
        res_mask = np.array([bool(lbl) for lbl in res_labels], dtype=bool)
        if res_mask.any():
            try:
                residues.ribbon_colors = _apply_rgb(
                    res_labels, res_mask, chain_colors, residues.ribbon_colors
                )
                residues.ring_colors = _apply_rgb(
                    res_labels, res_mask, chain_colors, residues.ring_colors
                )
            except (ValueError, IndexError) as exc:
                session.logger.warning(
                    f"[label-asym] colorbylabel: ribbon/ring color update "
                    f"failed ({exc}); atom colors applied but cartoons may "
                    f"not match"
                )

    n_atoms = int(atom_mask.sum())
    unique = len({lbl for lbl, m in zip(atom_labels, atom_mask) if m})
    session.logger.info(
        f"[label-asym] colorbylabel: colored {n_atoms} atoms across {unique} label chains"
    )


def register_commands(logger):
    """Register the ``colorbylabel`` command with ChimeraX's command registry."""
    from chimerax.core.commands import CmdDesc, ObjectsArg, register

    desc = CmdDesc(
        optional=[("objects", ObjectsArg)],
        synopsis="Color atoms by mmCIF label_asym_id (bychain analogue)",
    )
    register("colorbylabel", desc, colorbylabel, logger=logger)
