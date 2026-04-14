"""ChimeraX commands for the LabelAsym bundle.

Currently provides ``labelcolor``, a ``color bychain`` analogue that
groups residues by ``label_asym_id`` instead of ``auth_asym_id``. The
default palette matches ``chimerax.atomic.colors.chain_colors`` so heme
chains (label E/F/G/H in 4hhb) get visibly distinct colors from their
associated proteins (label A/B/C/D) while preserving the same palette
semantics users already know from ``color bychain``.
"""

from __future__ import annotations

try:
    from .hook import ATTR_NAME
except ImportError:  # tests import bare modules via sys.path
    from hook import ATTR_NAME


def _labels_of(items) -> list[str]:
    """Read ``label_asym_id`` per item, using "" when the attribute is unset."""
    return [getattr(it, ATTR_NAME, None) or "" for it in items]


def _apply_rgb(labels, mask, palette_fn, current_rgba):
    """Return a new RGBA array where ``mask`` entries are recolored by palette.

    Preserves the alpha channel of ``current_rgba`` and leaves masked-out rows
    untouched. ``palette_fn`` is called with the list of labels kept by
    ``mask`` and must return an ``(n, 4)`` uint8 array.
    """
    import numpy as np

    keep = [lbl for lbl, m in zip(labels, mask) if m]
    palette = palette_fn(keep)
    out = current_rgba.copy()
    idx = np.where(mask)[0]
    out[idx, :3] = palette[:, :3]
    return out


def labelcolor(session, objects=None):
    """Color atoms, ribbons, and rings by ``label_asym_id``."""
    import numpy as np

    from chimerax.atomic.colors import chain_colors
    from chimerax.core.objects import all_objects

    if objects is None:
        objects = all_objects(session)

    atoms = objects.atoms
    if len(atoms) == 0:
        session.logger.warning("[label-asym] labelcolor: no atoms selected")
        return

    atom_labels = _labels_of(atoms.residues)
    atom_mask = np.array([bool(lbl) for lbl in atom_labels], dtype=bool)
    if not atom_mask.any():
        session.logger.warning(
            "[label-asym] labelcolor: no residues have label_asym_id; "
            "open an mmCIF structure first"
        )
        return

    atoms.colors = _apply_rgb(atom_labels, atom_mask, chain_colors, atoms.colors)

    residues = objects.residues
    if len(residues) > 0:
        res_labels = _labels_of(residues)
        res_mask = np.array([bool(lbl) for lbl in res_labels], dtype=bool)
        if res_mask.any():
            residues.ribbon_colors = _apply_rgb(
                res_labels, res_mask, chain_colors, residues.ribbon_colors
            )
            residues.ring_colors = _apply_rgb(
                res_labels, res_mask, chain_colors, residues.ring_colors
            )

    n_atoms = int(atom_mask.sum())
    unique = len({lbl for lbl, m in zip(atom_labels, atom_mask) if m})
    session.logger.info(
        f"[label-asym] labelcolor: colored {n_atoms} atoms across {unique} label chains"
    )


def register_commands(logger):
    from chimerax.core.commands import CmdDesc, ObjectsArg, register

    desc = CmdDesc(
        optional=[("objects", ObjectsArg)],
        synopsis="Color atoms by mmCIF label_asym_id (bychain analogue)",
    )
    register("labelcolor", desc, labelcolor, logger=logger)
