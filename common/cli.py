"""Schema-driven CLI parsing for omnibenchmark module entrypoints (Python).

Reserved, overwrite-on-update path (``src/common/python/``) — see AGENTS.md.
The CLI an entrypoint accepts is defined as **data**, once, in
``src/common/schema/<interface>.json``, and built into a parser identically here
and in ``src/common/r/cli.R`` so Python and R entrypoints share one contract.

An *interface* is a named, versioned CLI contract owned by a benchmark; a module
"satisfies" it by carrying its schema and parsing against it. Schema shape::

    {
      "interface": "embedding",
      "version": "0.1.0",
      "benchmark": "omni-scrna/split-stages-plan",
      "args": [
        {"flag": "--name", "type": "string", "help": "...", "dest": "<optional>"},
        {"flag": "--solver", "type": "string", "choices": ["arpack", "randomized"]},
        ...
      ]
    }

Conventions: every arg is **required** (a run is reproducible from its
invocation line); unknown flags are rejected; ``dest`` defaults to the flag with
dots/dashes turned into ``_`` (``--pcas.tsv`` -> ``pcas_tsv``) unless the schema
overrides it. Types: ``path | string | integer | number``. An optional
``choices`` list restricts accepted values (an enum), like ``argparse``.

Layering — an interface is composed from up to three files in ``schema/``, each
contributing args, later layers winning per ``flag``:

  * ``_base.json``               — universal args every module gets
                                   (``--output_dir``, ``--name``); not a stage.
  * ``<interface>.json``         — the stage's I/O contract (benchmark-owned).
  * ``<interface>.extends.json`` — module-local extras/overrides (author-owned;
                                   never overwritten on update).

``parse_args("pca")`` discovers and merges all three by convention; a module that
carries only ``<interface>.json`` behaves exactly as before. Files starting with
``_`` or ending ``.extends.json`` are not stages (auto-pick skips them).

Usage in an entrypoint (in a rendered module the shared code is the ``common``
package)::

    from common.cli import parse_args
    args = parse_args("embedding")        # or parse_args() if the module has one schema
    # args.output_dir, args.name, args.pcas, args.clusters_truth
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

__version__ = "0.2.0"  # x-release-version — stamped from src/common/VERSION by `pixi run version`

def _find_schema_dir() -> Path:
    # Works under both layouts: rendered `common/cli.py` (schema is a sibling)
    # and the template's `common/python/cli.py` (schema is one level up).
    here = Path(__file__).resolve().parent
    for base in (here, here.parent):
        if (base / "schema").is_dir():
            return base / "schema"
    return here / "schema"


_SCHEMA_DIR = _find_schema_dir()
_TYPES = {"path": Path, "string": str, "integer": int, "number": float}


def common_version() -> str:
    """Version of the src/common shared code, so a module can report which copy
    of the boilerplate scaffolding it carries. Stamped from src/common/VERSION
    (single source of truth) — bump VERSION and run `pixi run version`."""
    return __version__


def _default_dest(flag: str) -> str:
    return flag.lstrip("-").replace(".", "_").replace("-", "_")


_BASE_SCHEMA = "_base"            # universal args, merged first; not a stage
_EXTENDS_SUFFIX = ".extends.json"  # module-local overlay for <interface>


def _is_stage_schema(path: Path) -> bool:
    """A stage schema is a plain ``<interface>.json`` — not the universal base
    (``_*``) and not an overlay (``*.extends.json``)."""
    return not path.name.startswith("_") and not path.name.endswith(_EXTENDS_SUFFIX)


def _read_schema(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def _merge_args(parent: list[dict], child: list[dict]) -> list[dict]:
    """Overlay ``child`` onto ``parent``: a child arg with the same ``flag``
    overrides (keeping the parent's position); a new ``flag`` is appended."""
    merged = {a["flag"]: a for a in parent}
    for a in child:
        merged[a["flag"]] = a
    return list(merged.values())


def load_interface(interface: str | None = None, schema_dir: Path = _SCHEMA_DIR) -> dict:
    """Load and compose an interface: ``_base`` (universal) -> ``<interface>``
    (the stage contract) -> ``<interface>.extends`` (module-local overrides),
    later layers winning per ``flag``. If ``interface`` is None, auto-pick the
    sole *stage* schema (``_*`` and ``*.extends.json`` are not stages)."""
    if interface is None:
        stages = [p for p in sorted(schema_dir.glob("*.json")) if _is_stage_schema(p)]
        if len(stages) != 1:
            raise SystemExit(
                f"specify an interface; found {len(stages)} stage schemas in {schema_dir}")
        interface = stages[0].name[: -len(".json")]

    stage_path = schema_dir / f"{interface}.json"
    if not stage_path.is_file():
        raise SystemExit(f"interface schema not found: {stage_path}")
    spec = _read_schema(stage_path)

    args: list[dict] = []
    base_path = schema_dir / f"{_BASE_SCHEMA}.json"
    if base_path.is_file():
        args = _merge_args(args, _read_schema(base_path).get("args", []))
    args = _merge_args(args, spec.get("args", []))
    extends_path = schema_dir / f"{interface}{_EXTENDS_SUFFIX}"
    if extends_path.is_file():
        args = _merge_args(args, _read_schema(extends_path).get("args", []))

    spec["args"] = args
    return spec


def build_parser(schema: dict) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=f"{schema.get('interface', 'module')} interface "
                    f"v{schema.get('version', '?')}")
    for arg in schema["args"]:
        parser.add_argument(
            arg["flag"],
            required=True,
            dest=arg.get("dest", _default_dest(arg["flag"])),
            type=_TYPES.get(arg["type"], str),
            choices=arg.get("choices"),  # None -> unconstrained
            help=arg.get("help", ""),
        )
    return parser


def parse_args(interface: str | None = None, argv=None) -> argparse.Namespace:
    return build_parser(load_interface(interface)).parse_args(argv)


if __name__ == "__main__":  # smoke / live demo
    a = parse_args()
    for k, v in sorted(vars(a).items()):
        print(f"{k}={v}")
