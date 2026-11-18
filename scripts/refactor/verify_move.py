"""Check that every top-level statement of the original module occurs verbatim
in the new public module or in one module of the internal package.

Usage: python verify_move.py <original.py> <public.py> <package_dir> [old=>new ...]
"""
import ast
import glob
import os
import sys

orig, public, pkg = sys.argv[1:4]
rewrites = [a.split("=>") for a in sys.argv[4:]]
src = open(orig, "rb").read().decode("utf-8").replace("\r\n", "\n")
for old, new in rewrites:
    src = src.replace(old, new)
tree = ast.parse(src)
targets = {public: open(public, "rb").read().decode("utf-8").replace("\r\n", "\n")}
for path in glob.glob(os.path.join(pkg, "*.py")):
    targets[path] = open(path, "rb").read().decode("utf-8").replace("\r\n", "\n")
missing = 0
placed = {}
for node in tree.body:
    if isinstance(node, (ast.Import, ast.ImportFrom)):
        continue
    start = node.decorator_list[0].lineno if getattr(node, "decorator_list", None) else node.lineno
    text = "\n".join(src.splitlines()[start - 1:node.end_lineno])
    homes = [p for p, t in targets.items() if text in t]
    if len(homes) != 1:
        missing += 1
        print("NOT FOUND EXACTLY ONCE:", start, homes)
    else:
        placed[homes[0]] = placed.get(homes[0], 0) + 1
total_new = sum(len(ast.parse(t).body) for t in targets.values())
print("statements checked:", sum(placed.values()), "missing:", missing)
for p in sorted(placed):
    print(f"  {os.path.basename(p):18s} {placed[p]}")
sys.exit(1 if missing else 0)
