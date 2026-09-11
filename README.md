# Scientific-Figure Normalization System

A local, deterministic toolchain for auditing and normalizing scientific
figures authored in **Adobe Illustrator** (with MATLAB-generated panels).

This repository's first component is a **read-only Illustrator analyzer** that
exports a machine-readable manifest of an `.ai` document — every artboard,
layer, group, text run, path and image — so downstream tooling (style linting,
font/size normalization, MATLAB cross-checks) can reason about a figure without
opening Illustrator interactively.

> This is **not** an autonomous GUI agent. It is a single deterministic
> ExtendScript that inspects `app.activeDocument` and writes JSON. It never
> modifies your document.

## Contents

| File | Purpose |
|------|---------|
| `illustrator/analyze_document.jsx` | Read-only document analyzer (ExtendScript / JSX). |
| `schemas/figure_manifest.schema.json` | JSON Schema (draft-07) describing the manifest. |
| `tests/sample_expected_manifest.json` | Example manifest for a simple Arial/RGB figure. |
| `README.md` | This file. |

## What the analyzer exports

For each run it writes **two files** into a directory you choose:

- `<docname>_manifest.json` — full structured manifest (UTF-8).
- `<docname>_summary.txt` — human-readable object counts + font usage.

The manifest captures:

- **Document**: filename, full path, Illustrator version, color space, ruler
  units, and every artboard with its `artboardRect`, width and height.
- **Hierarchy**: layers, nested sublayers, groups, compound paths and clipping
  groups, preserving parent → child relationships.
- **Every page item**: `uuid` (Illustrator 24.0+), `typename`, `name`, `note`,
  locked/hidden/editable state, `geometricBounds` / `visibleBounds` /
  `controlBounds`, `position`, `width`, `height`, `rotation` (see caveat below)
  and a parent descriptor.
- **TextFrameItem**: `contents`, text-frame kind, orientation, per-paragraph
  justification, and character **runs** coalesced by identical formatting —
  each run carrying font name/family/style, size, baseline shift, horizontal
  and vertical scale, tracking, capitalization, and fill/stroke color.
- **PathItem**: closed/filled/stroked flags, fill and stroke color, stroke
  width, path-point count, and clipping state.
- **RasterItem / PlacedItem**: embedded vs linked, source file path (when
  available), and effective bounds.
- **Statistics**: object counts by type, max nesting depth, and font usage.
- **Warnings**: every property that could not be read is recorded rather than
  aborting the run.

## Running the analyzer in Illustrator

The script targets Illustrator via `#target illustrator` and works with any
recent version (uuid requires Illustrator 24.0 / 2020 or newer; older versions
still run, with `uuid` reported as `null`).

1. **Open your figure** in Adobe Illustrator and make it the active document.
   (The script analyzes `app.activeDocument`.)
2. Run the script using **either** method:

   **A. From the menu (recommended)**
   - `File` ▸ `Scripts` ▸ `Other Script…`
     (macOS: `File > Scripts > Other Script…`;
     Windows: `File > Scripts > Other Script…`)
   - Select `illustrator/analyze_document.jsx`.

   **B. Install it as a listed script**
   - Copy `analyze_document.jsx` into Illustrator's Scripts folder, then
     restart Illustrator. It appears under `File ▸ Scripts ▸ analyze_document`.
     - **Windows:** `C:\Program Files\Adobe\Adobe Illustrator <ver>\Presets\en_US\Scripts\`
     - **macOS:** `/Applications/Adobe Illustrator <ver>/Presets/en_US/Scripts/`

   **C. From an editor (development)**
   - Open the file in the **VS Code ExtendScript Debugger** extension (or the
     legacy ExtendScript Toolkit), set the target application to Illustrator,
     and run.

3. When prompted, **choose an output directory**. Nothing is written until you
   pick one; cancelling the dialog aborts with no changes.
4. If an output file already exists, you are **asked to confirm** before it is
   overwritten. Declining aborts without writing.
5. On success an alert reports item and warning counts and the two output
   paths.

### Command-line note (optional, macOS)

You can drive Illustrator headlessly on macOS via AppleScript:

```bash
osascript -e 'tell application "Adobe Illustrator" to do javascript (read (POSIX file "/absolute/path/to/illustrator/analyze_document.jsx") as «class utf8»)'
```

On Windows the equivalent is the `Adobe Illustrator` COM automation
`DoJavaScriptFile`. The interactive `File ▸ Scripts` route above is the
supported, portable path.

## Illustrator scripting constraints handled by the script

These are the object-model limitations the analyzer defends against — worth
knowing if you extend it:

- **No `JSON` object.** ExtendScript uses an ES3-era engine, so the script
  ships its own JSON encoder and string-escaper (correct `\uXXXX`, control,
  and line/paragraph-separator escaping). Only classic `for` loops are used;
  `Array.forEach` / `indexOf` / `String.trim` are avoided.
- **`PageItem.uuid` is version-gated** (Illustrator 24.0 / 2020+). It is read
  through a try/catch and reported as `null` on older versions.
- **No generic `rotation`.** The DOM does not expose a per-item rotation angle
  (rotation lives in an unreadable transformation matrix), so `rotation` is
  always `null` and a single warning is emitted. This is intentional, not a
  bug.
- **`controlBounds` and `editable`** are not present on every item/version;
  reads are guarded and fall back to `null`.
- **`container.pageItems`** is used for *direct* children; nested groups and
  compound paths are recursed explicitly to preserve true hierarchy.
- **Colors are polymorphic** (RGB / CMYK / Gray / Lab / Spot / Gradient /
  Pattern / NoColor). Each is serialized by its `typename` with only the fields
  valid for that subtype.
- **Character attributes are undefined until set.** Runs are built by reading
  each character's concrete attributes and coalescing identical neighbors,
  which avoids `(mixed)` sentinels from range-level reads.
- **Failures are localized.** Any single unreadable property records a warning
  and the analysis continues; the whole run never aborts on one bad value.

## Validating a manifest

The manifest conforms to `schemas/figure_manifest.schema.json` (JSON Schema
draft-07). The schema is intentionally permissive (`additionalProperties: true`)
so future script versions can add fields without breaking validation.

```bash
python -m pip install jsonschema
python - <<'PY'
import json, jsonschema
schema = json.load(open("schemas/figure_manifest.schema.json"))
inst   = json.load(open("tests/sample_expected_manifest.json"))
jsonschema.validate(inst, schema)
print("manifest is valid")
PY
```

`tests/sample_expected_manifest.json` shows the exact shape produced for a
simple Arial/RGB figure containing a title, a mixed-font label, an axis path,
and a raster inside a clipping group.
