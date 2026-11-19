"""Split a module into an internal package plus a re-exporting public module.

Usage (from the project root): python scripts/refactor/split_module.py <plan.py> [--check]

The plan module defines:
  SOURCE        path of the module to split (relative to the project root)
  PACKAGE       dotted name of the internal package, e.g. "analyses._blast"
  PACKAGE_DIR   directory of the internal package
  PACKAGE_DOC   docstring of the package __init__
  MODULES       ordered dict  name -> (docstring, [top-level names])
  EXTRA         dict  synthetic name -> (module, source text inserted at the top of that module)
  REWRITES      list of (old, new) text replacements applied to the source first
  FACADE_NOTE   comment block placed after the public module's docstring

Every top-level statement moves byte for byte (after REWRITES) to the module
its names are assigned to; comments and blank lines before a statement move
with it. Statements that define no name (``if __name__ == "__main__":``)
stay in the public module. Imports are generated from the names each module
uses, and the module graph is checked for cycles before anything is written.
"""

from __future__ import annotations

import ast
import importlib.util
import os
import sys
from collections import defaultdict

ROOT = os.getcwd()


def load_plan(path):
    spec = importlib.util.spec_from_file_location("plan", path)
    plan = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plan)
    return plan


def defined_names(node):
    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
        return [node.name]
    out = []
    if isinstance(node, ast.Assign):
        for t in node.targets:
            for x in ast.walk(t):
                if isinstance(x, ast.Name):
                    out.append(x.id)
    elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
        out.append(node.target.id)
    return out


def annotation_nodes(node):
    for x in ast.walk(node):
        if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef)):
            args = x.args
            for a in args.posonlyargs + args.args + args.kwonlyargs + [args.vararg, args.kwarg]:
                if a is not None and a.annotation is not None:
                    yield a.annotation
            if x.returns is not None:
                yield x.returns
        elif isinstance(x, ast.AnnAssign):
            yield x.annotation


def used_names(node):
    names = {x.id for x in ast.walk(node) if isinstance(x, ast.Name) and not isinstance(x.ctx, ast.Store)}
    # Names a function imports itself are not module-level uses (Phase 9:
    # orf's finder imports DNA_IUPAC from core.alphabets inside a function).
    if not isinstance(node, (ast.Import, ast.ImportFrom)):
        for inner in ast.walk(node):
            if isinstance(inner, (ast.Import, ast.ImportFrom)) and inner is not node:
                names -= {(alias.asname or alias.name).split(".")[0] for alias in inner.names}
    for ann in annotation_nodes(node):
        for x in ast.walk(ann):
            if isinstance(x, ast.Constant) and isinstance(x.value, str):
                try:
                    expr = ast.parse(x.value, mode="eval")
                except SyntaxError:
                    continue
                names |= {y.id for y in ast.walk(expr) if isinstance(y, ast.Name)}
    return names


def main():
    plan = load_plan(sys.argv[1])
    check_only = "--check" in sys.argv
    src_path = os.path.join(ROOT, plan.SOURCE)
    raw = open(src_path, "rb").read().decode("utf-8")
    eol = "\r\n" if "\r\n" in raw else "\n"
    for old, new in plan.REWRITES:
        count = raw.count(old)
        if count == 0:
            sys.exit(f"rewrite not found: {old!r}")
        print(f"rewrite x{count}: {old!r} -> {new!r}")
        raw = raw.replace(old, new)
    lines = raw.splitlines(keepends=True)
    tree = ast.parse(raw)

    # Header: docstring and imports up to the first other statement.
    body = list(tree.body)
    doc_node = body[0] if isinstance(body[0], ast.Expr) and isinstance(getattr(body[0], "value", None), ast.Constant) else None
    imports = {}  # local name -> import line text
    header_end = doc_node.end_lineno if doc_node else 0
    rest = []
    for node in body[1 if doc_node else 0:]:
        if isinstance(node, (ast.Import, ast.ImportFrom)) and not rest:
            header_end = node.end_lineno
            if isinstance(node, ast.ImportFrom) and node.module == "__future__":
                continue
            for a in node.names:
                local = a.asname or a.name.split(".")[0]
                if isinstance(node, ast.ImportFrom):
                    imports[local] = ("from", node.module, a.name, a.asname)
                else:
                    imports[local] = ("import", a.name, None, a.asname)
        else:
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                sys.exit(f"import after code at line {node.lineno}; handle it in the plan")
            rest.append(node)

    assign = {}
    for mod, (_doc, names) in plan.MODULES.items():
        for n in names:
            if n in assign:
                sys.exit(f"{n} assigned twice")
            assign[n] = mod
    for n, (mod, _text) in plan.EXTRA.items():
        assign[n] = mod

    segments = defaultdict(list)   # module -> list of (text, node)
    facade_segments = []
    prev_end = header_end
    top_names = set(plan.EXTRA)
    for node in rest:
        start = prev_end + 1
        text = "".join(lines[start - 1:node.end_lineno])
        prev_end = node.end_lineno
        names = defined_names(node)
        top_names.update(names)
        if not names:
            is_main = isinstance(node, ast.If) and "__name__" in ast.unparse(node.test)
            if is_main:
                facade_segments.append((text, node))
                continue
            # A statement that only acts on existing names (such as
            # ``X.update(...)``) belongs with the module that defines them.
            owners = {assign[n] for n in used_names(node) if n in assign}
            if len(owners) != 1:
                sys.exit(f"cannot place statement at line {node.lineno}: owners {owners}")
            owner = owners.pop()
            print(f"placed non-defining statement at line {node.lineno} in {owner}")
            segments[owner].append((text, node))
            continue
        if set(names) <= set(getattr(plan, "FACADE_NAMES", ())):
            facade_segments.append((text, node))
            continue
        mods = {assign.get(n) for n in names}
        if None in mods:
            sys.exit(f"unassigned: {[n for n in names if n not in assign]} (line {node.lineno})")
        if len(mods) != 1:
            sys.exit(f"statement at line {node.lineno} defines names for several modules: {names}")
        segments[mods.pop()].append((text, node))
    tail = "".join(lines[prev_end:])
    if tail.strip() and not all(line.strip().startswith("#") or not line.strip() for line in tail.splitlines()):
        sys.exit("trailing code after the last statement")
    facade_tail = tail.strip("\r\n")
    missing = set(assign) - top_names
    if missing:
        sys.exit(f"assigned but not defined: {sorted(missing)}")

    # Dependencies.
    needs = {}
    graph = defaultdict(set)
    for mod in plan.MODULES:
        used = set()
        defined = set(plan_extra_names(plan, mod))
        for _text, node in segments[mod]:
            used |= used_names(node)
            defined |= set(defined_names(node))
        if mod in {m for m, _ in plan.EXTRA.values()}:
            for n, (m, text) in plan.EXTRA.items():
                if m == mod:
                    used |= used_names(ast.parse(text))
        internal = {n for n in used if n in assign and n not in defined}
        external = {n for n in used if n in imports and n not in defined and n not in assign}
        needs[mod] = (internal, external)
        for n in internal:
            graph[mod].add(assign[n])
    cycles = find_cycles(graph)
    for mod in plan.MODULES:
        print(f"{mod:14s} <- {sorted(graph[mod])}  ({len(segments[mod])} statements)")
    if cycles:
        sys.exit("CYCLES: " + "; ".join(" -> ".join(c) for c in cycles))
    if check_only:
        return

    os.makedirs(os.path.join(ROOT, plan.PACKAGE_DIR), exist_ok=True)
    write(os.path.join(ROOT, plan.PACKAGE_DIR, "__init__.py"), f'"""{plan.PACKAGE_DOC}"""\n', "\n")  # new files: LF
    for mod, (doc, _names) in plan.MODULES.items():
        internal, external = needs[mod]
        out = [f'"""{doc}"""\n', "\n", "from __future__ import annotations\n", "\n"]
        out += import_lines(external, imports)
        if external and internal:
            out.append("\n")
        by_mod = defaultdict(list)
        for n in internal:
            by_mod[assign[n]].append(n)
        for m in plan.MODULES:
            if m in by_mod:
                out += from_lines(f"{plan.PACKAGE}.{m}", sorted(by_mod[m], key=sort_key))
        for n, (m, text) in plan.EXTRA.items():
            if m == mod:
                out += ["\n", text if text.endswith("\n") else text + "\n"]
        body_text = "".join(t for t, _ in segments[mod]).lstrip("\r\n")
        out += ["\n", "\n", body_text]
        write(os.path.join(ROOT, plan.PACKAGE_DIR, f"{mod}.py"), "".join(out), "\n")

    # Public module.
    doc_text = "".join(lines[:doc_node.end_lineno]) if doc_node else ""
    out = [doc_text.replace("\r\n", "\n"), "\n", plan.FACADE_NOTE, "\n", "from __future__ import annotations\n", "\n"]
    facade_used = set()
    for _t, node in facade_segments:
        facade_used |= used_names(node)
    out += import_lines({n for n in facade_used if n in imports and n not in assign}, imports)
    for mod, (_doc, names) in plan.MODULES.items():
        exported = [n for n in names] + [n for n, (m, _t) in plan.EXTRA.items() if m == mod and n in plan.EXPORT_EXTRA]
        out += from_lines(f"{plan.PACKAGE}.{mod}", sorted(exported, key=sort_key))
    for text, _node in facade_segments:
        out += ["\n", "\n", text.lstrip("\r\n").replace("\r\n", "\n")]
    if facade_tail:
        out += ["\n", "\n", facade_tail.replace("\r\n", "\n"), "\n"]
    write(os.path.join(ROOT, plan.SOURCE), "".join(out), eol)


def plan_extra_names(plan, mod):
    return [n for n, (m, _t) in plan.EXTRA.items() if m == mod]


def sort_key(name):
    return (name.lstrip("_").lower(), name)


def import_lines(names, imports):
    plain = sorted({imports[n] for n in names if imports[n][0] == "import"})
    froms = defaultdict(list)
    for n in names:
        kind, module, name, asname = imports[n]
        if kind == "from":
            froms[module].append(name if not asname else f"{name} as {asname}")
    out = []
    for _k, module, _n, asname in plain:
        out.append(f"import {module}" + (f" as {asname}" if asname else "") + "\n")
    for module in sorted(froms):
        out += from_lines(module, sorted(froms[module], key=sort_key))
    return out


def from_lines(module, names):
    line = f"from {module} import {', '.join(names)}\n"
    if len(line) <= 100:
        return [line]
    return [f"from {module} import (\n"] + [f"    {n},\n" for n in names] + [")\n"]


def write(path, text, eol):
    text = text.replace("\r\n", "\n").replace("\n", eol)
    with open(path, "wb") as handle:
        handle.write(text.encode("utf-8"))
    print("wrote", os.path.relpath(path, ROOT))


def find_cycles(graph):
    cycles, color = [], {}

    def dfs(u, stack):
        color[u] = 1
        for v in sorted(graph[u]):
            if color.get(v) == 1:
                cycles.append(stack + [u, v])
            elif not color.get(v):
                dfs(v, stack + [u])
        color[u] = 2

    for u in list(graph):
        if not color.get(u):
            dfs(u, [])
    return cycles


main()
