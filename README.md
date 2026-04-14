# ChimeraX-LabelAsym

ChimeraX bundle that exposes the mmCIF `label_asym_id` as a residue attribute,
so you can select residues by the PDBx-canonical chain ID in addition to the
author-assigned chain ID (`auth_asym_id`) that ChimeraX uses by default.

## Why

In mmCIF:

- `auth_asym_id` — author-assigned chain (what `/A` in ChimeraX selects).
- `label_asym_id` — PDB-canonical chain.

They can differ. For example, in `4hhb` the heme of chain A has
`auth_asym_id=A` but `label_asym_id=E`. With only `auth_asym_id` there is no
way to pick exactly that heme (and nothing else) through `/A`.

## Usage

After installing the bundle, any structure opened from mmCIF metadata gets a
`label_asym_id` attribute on every residue automatically:

```
open 4hhb format mmcif
select ::label_asym_id="E"   # just the heme of chain A
color  ::label_asym_id="A" red
```

The built-in `/A` syntax continues to work as before (author chain).

## Build & Install

```bash
echi build
echi install
```

## Test

```bash
echi run --script scripts/smoke.cxc
```

## Future work

- Short syntax such as `/A'` that means "label_asym_id=A" — requires extending
  ChimeraX's atom-spec parser.

## License

MIT
