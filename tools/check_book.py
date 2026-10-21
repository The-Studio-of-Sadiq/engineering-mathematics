"""Consistency checks for the whole book. Exit code 1 if anything fails.

    python tools/check_book.py          # report problems only
    python tools/check_book.py -v       # also list every file with its word count

Checks every Markdown file (outside node_modules/.git) for:
  - unbalanced $$ display math and ``` code fences
  - stray control characters (shell heredocs turn \\f, \\b, \\t in LaTeX into them)
  - the "^*{" corruption artifact
  - non-portable LaTeX (\\tag, \\label, \\eqref, \\newcommand, align) outside 00-front/
  - exercise/answer headings
  - relative .md and image links that point at missing files
  - chapters: section headings numbered for another chapter, or duplicated
  - track chapters (26-99): the seven-part arc, at least one figure, at least two [!mystery] callouts
"""

import collections
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKIP = {"node_modules", ".git", ".obsidian"}
TRACK_SECTIONS = ["The hook", "The big idea", "The mathematics", "The algorithm",
                  "The experiment", "The challenge", "Loose threads"]


def problems_in(path):
    s = path.read_text(encoding="utf-8")
    rel = path.relative_to(ROOT).as_posix()
    nocode = re.sub(r"```.*?```", "", s, flags=re.S)
    nocode = re.sub(r"``[^`\n]+``|`[^`\n]+`", "", nocode)  # inline code too: it holds examples, not content
    out = []
    if nocode.count("$$") % 2:
        out.append("odd number of $$")
    if s.count("```") % 2:
        out.append("odd number of ``` fences")
    ctrl = sorted({hex(ord(c)) for c in s if ord(c) < 32 and c not in "\n\r\t"})
    if ctrl:
        out.append("control characters " + ",".join(ctrl))
    if "^*{" in s:
        out.append(f"'^*{{' artifact x{s.count('^*{')}")
    if not rel.startswith("00-front/"):
        for pat, name in [(r"\\tag\{", r"\tag"), (r"\\label\{", r"\label"), (r"\\eqref", r"\eqref"),
                          (r"\\newcommand", r"\newcommand"), (r"\\begin\{align\*?\}", "align")]:
            if re.search(pat, nocode):
                out.append(f"uses {name}")
    if re.search(r"^#+ .*\b(Exercises?|Answers)\b", s, re.M):
        out.append("exercise/answers heading")
    for link in re.findall(r"\]\(([^)\s]+)\)", nocode):
        target = link.split("#")[0]
        if not target or re.match(r"^[a-z]+:", target) or not target.endswith((".md", ".svg", ".png", ".pdf")):
            continue
        if not (path.parent / target).resolve().exists():
            out.append(f"broken link {link}")
    m = re.match(r"(\d\d|[A-Z])-", path.name)
    if m and path.parent.name != "00-front":
        key = m.group(1).lstrip("0")
        heads = re.findall(r"^## (\d+|[A-Z])\.(\d+)", nocode, re.M)
        foreign = sorted({f"{a}.{b}" for a, b in heads if a != key})
        if foreign:
            out.append("sections from another chapter: " + ", ".join(foreign[:5]))
        dup = [f"{a}.{b}" for (a, b), n in collections.Counter(heads).items() if n > 1]
        if dup:
            out.append("duplicate sections: " + ", ".join(dup))
        if rel.startswith("tracks/"):
            titles = re.findall(r"^## \d+\.\d+ (.+)$", nocode, re.M)
            missing = [t for t in TRACK_SECTIONS if not any(x.strip().startswith(t) for x in titles)]
            if missing:
                out.append("missing track sections: " + ", ".join(missing))
            if not re.search(r"!\[[^\]]*\]\([^)]*images/", s):
                out.append("no figure")
            if s.count("[!mystery]") < 2:
                out.append("fewer than two [!mystery] callouts")
    return out, len(s.split())


def main():
    verbose = "-v" in sys.argv
    files = sorted(p for p in ROOT.rglob("*.md") if not SKIP.intersection(p.parts))
    bad, words = 0, 0
    for p in files:
        probs, w = problems_in(p)
        words += w
        rel = p.relative_to(ROOT).as_posix()
        if probs:
            bad += 1
            print(f"{rel}: " + "; ".join(probs))
        elif verbose:
            print(f"{rel}: ok ({w} words)")
    print(f"{len(files)} files, {words:,} words, {bad} with problems")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
