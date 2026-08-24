#!/usr/bin/env python3
"""Pull stage-input fixtures for THIS module out of the benchmark's latest fixture run.

Three steps (per request):
  1. download the `fixture-run-out` artifact from the latest successful fixture
     run of the master plan (the benchmark named in omnibenchmark.yaml);
  2. identify which stages this module implements — plan stages whose module
     points at this repo AND uses one of our omnibenchmark.yaml entrypoints
     (the plan IS the "implements" annotation we don't have locally);
  3. for each such stage, resolve its declared inputs to filenames and copy one
     sample of each out of the run tree into <dest>/<stage>/.

The plan is the source of truth for input->file wiring: the vendored
src/common/schema/*.json carry only the input *flags*, not the filenames.

    pixi run python ci/extract-fixtures.py --stage PCA
    pixi run python ci/extract-fixtures.py --selftest   # logic check, no network

ponytail: downloads the whole 1.5 GB artifact to grab a few small files — the
only artifact that exists. Switch to the contract-fixtures Release bundle once
the plan publishes it (small, stable); see docs/stage-contract-test-data.md.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml

# omnibenchmark bookkeeping dirs in a run tree; never a real output.
META = {".snakemake", ".modules", ".metadata", ".envs", ".logs", ".git"}


def sh(*cmd: str) -> str:
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {' '.join(cmd)}\n{r.stderr.strip()}")
    return r.stdout.strip()


def repo_slug() -> str:
    """owner/name for this checkout (strips ssh-host aliases like github-traven)."""
    url = sh("git", "remote", "get-url", "origin")
    m = re.search(r"[:/]([^/:]+/[^/]+?)(?:\.git)?$", url)
    if not m:
        sys.exit(f"could not parse repo slug from remote: {url}")
    return m.group(1)


def fetch_plan(repo: str, path: str, ref: str) -> dict:
    """The benchmark plan YAML, read at the ref this module targets."""
    raw = sh("gh", "api", f"repos/{repo}/contents/{path}?ref={ref}", "--jq", ".content")
    import base64
    return yaml.safe_load(base64.b64decode(raw))


def implements(plan: dict, slug: str, entrypoints: set[str]) -> dict[str, dict]:
    """stage_id -> {entrypoint, inputs} for stages this repo provides a module for."""
    found = {}
    for stage in plan.get("stages", []):
        for mod in stage.get("modules", []) or []:
            url = (mod.get("repository") or {}).get("url", "")
            ep = (mod.get("repository") or {}).get("entrypoint") or mod.get("id", "")
            if slug in url and (not entrypoints or ep in entrypoints):
                found[stage["id"]] = {
                    "entrypoint": ep,
                    "inputs": list(stage.get("inputs", []) or []),
                }
                break
    return found


def output_paths(plan: dict) -> dict[str, str]:
    """output-id -> path template, across every stage."""
    return {
        o["id"]: o["path"]
        for stage in plan.get("stages", [])
        for o in (stage.get("outputs", []) or [])
    }


def dataset_name(plan: dict, override: str | None) -> str:
    if override:
        return override
    for stage in plan.get("stages", []):
        if stage["id"] == "DATA":
            for mod in stage.get("modules", []) or []:
                for params in mod.get("parameters", []) or []:
                    if params.get("dataset_name"):
                        return params["dataset_name"]
    return "be1"  # fixtures/1-data default


def latest_run_id(plan_repo: str) -> str:
    out = sh("gh", "run", "list", "-R", plan_repo, "--workflow", "fixture-run.yml",
             "--branch", "main", "--status", "success", "-L", "1",
             "--json", "databaseId", "--jq", ".[0].databaseId")
    if not out:
        sys.exit(f"no successful fixture run found on {plan_repo}@main")
    return out


def download_run_tree(plan_repo: str, run_id: str, workdir: Path) -> Path:
    """Download + untar the fixture-run artifact; return the out_ci/ root."""
    sh("gh", "run", "download", run_id, "-R", plan_repo,
       "-n", "fixture-run-out", "-D", str(workdir))
    tarball = workdir / "fixture-run-out.tar.gz"
    sh("tar", "xzf", str(tarball), "-C", str(workdir))
    root = workdir / "out_ci"
    if not root.is_dir():
        sys.exit(f"artifact did not contain out_ci/ (looked in {workdir})")
    return root


def find_in_tree(root: Path, filename: str) -> Path | None:
    """First file named `filename` in the run tree (any valid producer will do)."""
    hits = [p for p in sorted(root.rglob(filename)) if not set(p.parts) & META]
    return hits[0] if hits else None


def resolve_inputs(stage: dict, paths: dict[str, str], dataset: str) -> dict[str, str]:
    """input output-id -> concrete filename for this dataset."""
    resolved = {}
    for in_id in stage["inputs"]:
        tmpl = paths.get(in_id)
        if not tmpl:
            sys.exit(f"input '{in_id}' has no declared output path in the plan")
        resolved[in_id] = tmpl.replace("{dataset}", dataset)
    return resolved


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", help="stage id to extract (default: all this module implements)")
    ap.add_argument("--dataset", help="dataset name (default: from plan DATA stage)")
    ap.add_argument("--dest", type=Path, default=Path("fixtures"), help="output dir")
    ap.add_argument("--run-id", help="fixture run id (default: latest successful on plan main)")
    ap.add_argument("--from-dir", type=Path,
                    help="use an already-extracted out_ci/ instead of downloading")
    ap.add_argument("--selftest", action="store_true", help="run logic check, no network")
    a = ap.parse_args()

    if a.selftest:
        return _selftest()

    omni = yaml.safe_load(Path("omnibenchmark.yaml").read_text())
    bench = omni["benchmarks"][0]
    entrypoints = set(omni.get("entrypoints", {}))
    slug = repo_slug()

    plan = fetch_plan(_plan_slug(bench["repo"]), bench["plan"], bench["ref"])
    stages = implements(plan, slug, entrypoints)
    if not stages:
        sys.exit(f"{slug} implements no stages in {bench['repo']}")
    print(f"{slug} implements: {', '.join(sorted(stages))}")

    targets = [a.stage] if a.stage else list(stages)
    for t in targets:
        if t not in stages:
            sys.exit(f"stage '{t}' not implemented by this module (have: {', '.join(stages)})")

    paths = output_paths(plan)
    dataset = dataset_name(plan, a.dataset)

    if a.from_dir:
        root = a.from_dir
        tmp = None
    else:
        plan_repo = _plan_slug(bench["repo"])
        run_id = a.run_id or latest_run_id(plan_repo)
        print(f"downloading fixture-run-out from {plan_repo} run {run_id} ...")
        tmp = Path(tempfile.mkdtemp(prefix="fixture-run-"))
        root = download_run_tree(plan_repo, run_id, tmp)

    try:
        for t in targets:
            wanted = resolve_inputs(stages[t], paths, dataset)
            out = a.dest / t
            out.mkdir(parents=True, exist_ok=True)
            for in_id, fname in wanted.items():
                src = find_in_tree(root, fname)
                if not src:
                    sys.exit(f"[{t}] input '{in_id}' ({fname}) not found in run tree")
                dst = out / fname
                shutil.copy2(src, dst)
                print(f"[{t}] {in_id}: {src.relative_to(root)} -> {dst}")
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)


def _plan_slug(repo_url: str) -> str:
    return re.sub(r"^https?://github\.com/|\.git$", "", repo_url)


def _selftest() -> None:
    plan = {
        "stages": [
            {"id": "DATA",
             "outputs": [{"id": "rawdata_h5ad", "path": "{dataset}.h5ad"}],
             "modules": [{"id": "data",
                          "repository": {"url": "https://github.com/omni-scrna/1-data"},
                          "parameters": [{"dataset_name": "be1"}]}]},
            {"id": "FEAT",
             "outputs": [{"id": "normalized_selected_h5", "path": "{dataset}_normalized_selected.h5"}],
             "modules": [{"id": "fe-scanpy",
                          "repository": {"url": "https://github.com/omni-scrna/scanpy",
                                         "entrypoint": "feat-select"}}]},
            {"id": "PCA",
             "inputs": ["normalized_selected_h5"],
             "outputs": [{"id": "pcas_tsv", "path": "{dataset}_pcas.tsv"}],
             "modules": [{"id": "pc-scanpy",
                          "repository": {"url": "https://github.com/omni-scrna/scanpy",
                                         "entrypoint": "pca"}},
                         {"id": "pc-scrapper",
                          "repository": {"url": "https://github.com/omni-scrna/scrapper",
                                         "entrypoint": "pca"}}]},
        ]
    }
    eps = {"pca", "feat-select", "knn"}
    impl = implements(plan, "omni-scrna/scanpy", eps)
    assert set(impl) == {"FEAT", "PCA"}, impl
    assert impl["PCA"]["entrypoint"] == "pca", impl
    # scrapper's pca module must not make us claim PCA via the wrong repo path
    assert "omni-scrna/scanpy" not in str(impl["PCA"].get("url", "")), impl

    assert dataset_name(plan, None) == "be1"
    paths = output_paths(plan)
    assert paths["normalized_selected_h5"] == "{dataset}_normalized_selected.h5"
    resolved = resolve_inputs(impl["PCA"], paths, "be1")
    assert resolved == {"normalized_selected_h5": "be1_normalized_selected.h5"}, resolved

    # repo slug parsing handles ssh host aliases
    assert re.search(r"[:/]([^/:]+/[^/]+?)(?:\.git)?$",
                     "git@github-traven.github.com:omni-scrna/scanpy.git").group(1) == "omni-scrna/scanpy"
    print("selftest OK")


if __name__ == "__main__":
    main()
