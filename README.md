# ChimeraX-LabelAsym

ChimeraX bundle that exposes the mmCIF `label_asym_id` as a residue attribute,
so you can select residues by the PDBx-canonical chain ID in addition to the
author-assigned chain ID (`auth_asym_id`) that ChimeraX uses by default.

## Why

In mmCIF:

- `auth_asym_id` — author-assigned chain (what `/A` in ChimeraX selects).
- `label_asym_id` — PDB-canonical chain.

They can differ. For example, in `4hhb` the four heme groups all share an
`auth_asym_id` with their host protein chain, but each has its own
`label_asym_id`:

| entity           | `auth_asym_id` | `label_asym_id` |
|------------------|:--------------:|:---------------:|
| α1 protein       | A              | A               |
| β1 protein       | B              | B               |
| α2 protein       | C              | C               |
| β2 protein       | D              | D               |
| heme of α1       | A              | E               |
| heme of β1       | B              | F               |
| heme of α2       | C              | G               |
| heme of β2       | D              | H               |

Without `label_asym_id` there is no way to pick exactly one heme through the
default `/A` syntax — it always carries the protein with it.

## Usage

After installing the bundle, any structure opened from mmCIF metadata gets a
`label_asym_id` attribute on every residue automatically:

```
open 4hhb format mmcif
# log: [label-asym] 4hhb: assigned label_asym_id to 801/801 residues
```

The bundle provides three ways to work with label chains:

### 1. Attribute selector (any label value)

```
color ::label_asym_id="A" dodger blue    # just the α1 protein
color ::label_asym_id="E" orange         # just its heme
view  ::label_asym_id="E"                # zoom to that heme
```

<p align="center">
  <img src="docs/04_label_E_zoom.png" width="520" alt="label_asym_id E = heme of chain A">
</p>

### 2. Short `la_<label>` selectors (registered on open)

Every unique `label_asym_id` assigned to at least one residue in a newly
opened structure gets a matching `la_<label>` selector (subject to the
character restrictions noted under Limitations), so you can skip the
attribute syntax:

```
select la_E              # heme of α1 (same as ::label_asym_id="E")
color  la_A dodger blue  # α1 protein
hide   la_I target a     # hide the first waters, etc.
```

Selectors are session-global: the same `la_E` matches residues across every
open model carrying that label. Opening a second structure that shares
labels reuses the existing selectors instead of duplicating them.

### 3. `labelcolor` — `color bychain` for label chains

Colors atoms, cartoons, and rings by `label_asym_id` using the same palette as
`color bychain`, so hemes (label E/F/G/H) get a colour distinct from their
associated proteins (label A/B/C/D):

```
labelcolor            # color everything by label_asym_id
labelcolor #1/A       # restrict to one author chain
```

<p align="center">
  <img src="docs/06_labelcolor_4hhb.png" width="520" alt="labelcolor 4hhb: protein and heme chains get distinct colors">
</p>

The built-in `/A` syntax continues to work as before (author chain).

### Every entity in its own colour (manual)

```
color ::label_asym_id="A" dodger blue
color ::label_asym_id="B" forest green
color ::label_asym_id="C" hot pink
color ::label_asym_id="D" gold
color ::label_asym_id="E" orange
color ::label_asym_id="F" red
color ::label_asym_id="G" purple
color ::label_asym_id="H" cyan
```

<p align="center">
  <img src="docs/05_all_labels.png" width="520" alt="all label_asym_id groups coloured independently">
</p>

Note: quotes around single-letter values are required in the attribute form
(e.g. `::label_asym_id="A"`). Without the quotes, ChimeraX's atom-spec parser
reads single letters such as `A`, `C`, or `G` as residue-code tokens and
rejects the selection. The `la_<label>` short selectors and `labelcolor`
command have no such quoting requirement.

## Build & Install

### Option A — ChimeraX's built-in `devel` / `toolshed` commands

Build a wheel directly from this source tree (ChimeraX does the packaging):

```bash
# From the repo root. Writes the wheel to ./dist/
/Applications/ChimeraX-*.app/Contents/MacOS/ChimeraX --nogui --exit \
    --cmd "devel build ."
```

Then install the wheel into your ChimeraX profile:

```bash
/Applications/ChimeraX-*.app/Contents/MacOS/ChimeraX --nogui --exit \
    --cmd "toolshed install dist/chimerax_labelasym-*.whl"
```

Or do both in one step from an already-running ChimeraX (Tools → Log,
or the command line in the GUI):

```
devel install /path/to/ChimeraX-LabelAsym
```

Restart ChimeraX afterward — `custom-init` runs at startup.

### Option B — [echidna](https://github.com/N283T/echidna)

[Echidna](https://github.com/N283T/echidna) is a Rust CLI that wraps the
above `devel` / `toolshed` invocations (auto-detects `CHIMERAX_PATH`, builds
and installs in one command, launches ChimeraX with the bundle):

```bash
echi build
echi install
# or all at once:
echi run
```

## Test

Unit tests (pure Python, no ChimeraX runtime required):

```bash
uv run --with pytest --no-project pytest
```

Integration smoke test (requires ChimeraX):

```bash
echi run --script scripts/smoke.cxc
```

## How it works

1. `custom-init = true` in `pyproject.toml` runs `initialize()` on ChimeraX
   startup.
2. `hook.install()` registers `Residue.label_asym_id` as a persistent string
   attribute and hooks the `ADD_MODELS` trigger.
3. On every new `AtomicStructure`, the hook pulls the `atom_site` CIF table
   (via `get_mmcif_tables_from_metadata`, falling back to re-reading the file
   with `get_cif_tables`) and builds a
   `(auth_asym_id, auth_seq_id, ins_code) → label_asym_id` map.
4. Each residue gets its `label_asym_id` attribute set. The key lookup is
   dual-scheme — a parallel `(label_asym_id, label_seq_id, "")` entry is
   stored alongside the auth entry so the same code works for structures
   opened with `prefer_auth=false`.
5. For every unique label assigned to at least one residue, a `la_<label>`
   selector is registered via `chimerax.core.commands.register_selector`.
   Registration is deduped across structures and skipped for labels that
   would violate ChimeraX's selector-name rules. `uninstall()` deregisters
   every selector this bundle created. Non-mmCIF structures are silently
   skipped.

## Limitations

- Only mmCIF-derived structures are annotated. PDB-format input does not carry
  `label_asym_id`.
- `la_<label>` selectors require labels composed of alphanumerics, `-`, `+`,
  or `_`. Unusual labels containing other punctuation are still reachable via
  the attribute form (`::label_asym_id="..."`).
- No `/A'` grammar extension — ChimeraX's atom-spec parser is not extensible
  from a bundle, so attribute selectors and `la_<label>` short names are used
  instead.

## License

MIT
