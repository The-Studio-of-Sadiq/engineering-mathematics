"""Repair Desmos-block syntax the obsidian-desmos plugin can't parse.

The plugin splits each equation line on "|" to read style options, so absolute values written
\\left|x\\right| are cut in half. Inside ```desmos-graph blocks only, this rewrites them as
\\operatorname{abs}\\left(x\\right), and turns double inequalities a<f<b (unsupported by Desmos
unless f is x or y) into f>a|f<b (an inequality plus a restriction), and makes products with
list literals explicit (\\frac{2\\pi}{7}[0,1,2] would otherwise index a number).

    python tools/fix_desmos.py            # fix in place, report counts
    python tools/fix_desmos.py --check    # report only; exit 1 if anything needs fixing
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BLOCK = re.compile(r"(```desmos-graph\n)(.*?)(```)", re.S)
ABS = re.compile(r"\\left\|(.*?)\\right\|")
NUM = r"-?\d+(?:\.\d+)?"
DOUBLE = re.compile(rf"^({NUM})\s*(<=?|\\le)\s*(.+?)\s*(<=?|\\le)\s*({NUM})$")


def fix_line(line):
    # \frac{a}{b}[0,1,2] reads as indexing in Desmos; make the list product explicit.
    line = re.sub(r"\}\[", r"}\\cdot[", line)
    head, sep, rest = line.partition("|") if not ABS.search(line) else (line, "", "")
    if ABS.search(line):
        while ABS.search(line):
            line = ABS.sub(lambda m: r"\operatorname{abs}\left(" + m.group(1) + r"\right)", line)
        head, sep, rest = line.partition("|")
    m = DOUBLE.match(head.strip())
    if m and m.group(3).strip() not in ("x", "y"):
        lo, op1, f, op2, hi = m.groups()
        flip = {"<": ">", "<=": ">=", "\\le": "\\ge"}
        head = f"{f}{flip[op1]}{lo}|{f}{op2}{hi}"
    return head + (sep + rest if sep else "")


def fix_block(body):
    lines = body.split("\n")
    out, in_eq = [], "---" not in [l.strip() for l in lines]
    for l in lines:
        if l.strip().startswith("---"):
            in_eq = True
            out.append(l)
            continue
        out.append(fix_line(l) if in_eq and l.strip() else l)
    return "\n".join(out)


def main():
    check = "--check" in sys.argv
    prefixes = tuple(a for a in sys.argv[1:] if not a.startswith("--"))  # optional path prefixes
    changed = 0
    for p in sorted(ROOT.rglob("*.md")):
        if any(part in ("node_modules", ".git") for part in p.parts):
            continue
        if prefixes and not p.relative_to(ROOT).as_posix().startswith(prefixes):
            continue
        text = p.read_text(encoding="utf-8")
        new = BLOCK.sub(lambda m: m.group(1) + fix_block(m.group(2)) + m.group(3), text)
        if new != text:
            n = sum(a != b for a, b in zip(text.splitlines(), new.splitlines()))
            changed += n
            print(f"{p.relative_to(ROOT).as_posix()}: {n} lines")
            if not check:
                p.write_text(new, encoding="utf-8")
    print(f"{changed} lines {'need fixing' if check else 'fixed'}")
    sys.exit(1 if check and changed else 0)


if __name__ == "__main__":
    main()
