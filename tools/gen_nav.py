"""Generate navigation pages from tools/book.json.

Writes:
  00-front/table-of-contents.md   (all three levels + appendix)
  00-front/syllabus-tracks.md     (what each track chapter explores, and its challenge)
  00-front/learning-paths.md      (curated routes through the tracks)
  tracks/<dir>/README.md          (one index page per track)

and, with --plan, prints a JSON reading plan (path, prev, next per chapter) used when writing chapters.

    python tools/gen_nav.py
"""

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BOOK = json.loads((ROOT / "tools" / "book.json").read_text(encoding="utf-8"))

CORE_PARTS = [
    ("I — Mathematical Foundations", "part-1-foundations"),
    ("II — Calculus", "part-2-calculus"),
    ("III — Linear Algebra", "part-3-linear-algebra"),
    ("IV — Differential Equations", "part-4-differential-equations"),
    ("V — Transforms and Mathematical Analysis", "part-5-transforms"),
    ("VI — Probability, Statistics, and Stochastic Mathematics", "part-6-probability-statistics"),
    ("VII — Numerical and Computational Mathematics", "part-7-numerical-computational"),
]


def title_of(path):
    first = path.read_text(encoding="utf-8").splitlines()[0]
    return first.lstrip("# ").split("—", 1)[-1].strip()


def core_chapters():
    out = {}
    for part, d in CORE_PARTS:
        for p in sorted((ROOT / d).glob("[0-9][0-9]-*.md")):
            out[int(p.name[:2])] = {"path": p.relative_to(ROOT).as_posix(), "title": title_of(p), "part": part}
    return out


def appendices():
    return [{"path": p.relative_to(ROOT).as_posix(), "letter": p.name[0], "title": title_of(p)}
            for p in sorted((ROOT / "appendix").glob("[A-Z]-*.md"))]


def track_chapters():
    """Reading order of all track chapters with paths, prev and next (paths relative to ROOT)."""
    order = []
    for t in BOOK["tracks"]:
        for c in t["chapters"]:
            order.append({**c, "track": t["id"], "track_title": t["title"], "level": t["level"],
                          "path": f"{t['dir']}/{c['n']:02d}-{c['slug']}.md", "index": f"{t['dir']}/README.md"})
    for i, c in enumerate(order):
        first = i == 0 or order[i - 1]["track"] != c["track"]
        last = i == len(order) - 1 or order[i + 1]["track"] != c["track"]
        c["prev"] = c["index"] if first else order[i - 1]["path"]
        c["next"] = (order[i + 1]["index"] if i + 1 < len(order) else "00-front/table-of-contents.md") if last \
            else order[i + 1]["path"]
    return order


def rel(src, dst):
    return os.path.relpath(ROOT / dst, (ROOT / src).parent).replace(os.sep, "/")


def write(path, text):
    (ROOT / path).parent.mkdir(parents=True, exist_ok=True)
    (ROOT / path).write_text(text, encoding="utf-8")
    print("wrote", path)


def toc(core, chapters, apps):
    src = "00-front/table-of-contents.md"
    L = ["# Table of Contents", "",
         f"[← README]({rel(src, 'README.md')}) · [Notation & Conventions](notation-and-conventions.md) · "
         "[Learning Paths](learning-paths.md)", "",
         "The book has three levels. **Level 1** is the engineering-mathematics core. **Level 2** is nine "
         "specialization tracks, each with its own prerequisite map. **Level 3** is the research frontier. "
         "Nobody needs to read it front to back: pick a track or a [learning path](learning-paths.md).", "",
         f"![The Mathematics Universe]({rel(src, 'images/universe-map.svg')})", "",
         "## Level 1 — Core", ""]
    for part, d in CORE_PARTS:
        L.append(f"### Part {part}")
        L += [f"{n}. [{c['title']}]({rel(src, c['path'])})" for n, c in core.items() if c["path"].startswith(d)]
        L.append("")
    for level in (2, 3):
        name = next(l["name"] for l in BOOK["levels"] if l["level"] == level)
        L += [f"## Level {level} — {name}", ""]
        for t in [t for t in BOOK["tracks"] if t["level"] == level]:
            L.append(f"### [Track {t['id']} — {t['title']}]({rel(src, t['dir'] + '/README.md')})")
            L += [f"- Ch. {c['n']} — [{c['title']}]({rel(src, c['path'])})" for c in chapters if c["track"] == t["id"]]
            L.append("")
    L += ["## Appendix — Derivations and Proofs", ""]
    L += [f"- {a['letter']}. [{a['title']}]({rel(src, a['path'])})" for a in apps]
    return "\n".join(L) + "\n"


def track_index(t, chapters, core):
    src = f"{t['dir']}/README.md"
    level = next(l["name"] for l in BOOK["levels"] if l["level"] == t["level"])
    mine = [c for c in chapters if c["track"] == t["id"]]
    roadmap = rel(src, f"images/roadmap-track-{t['id'].lower()}.svg")
    L = [f"# Track {t['id']} — {t['title']}", "",
         f"[← Table of Contents]({rel(src, '00-front/table-of-contents.md')}) · "
         f"[Start: Ch. {mine[0]['n']} →]({rel(src, mine[0]['path'])})", "",
         f"> **Level {t['level']} — {level}.** {t['blurb']}", "",
         f"![Track {t['id']} roadmap]({roadmap})", "",
         f"*Suggested reading order for Track {t['id']}. The grey stations on the left are the core chapters it leans on.*", "",
         "## Before you start", "",
         "This track assumes you are comfortable with these core chapters:", ""]
    L += [f"- Ch. {n} — [{core[n]['title']}]({rel(src, core[n]['path'])})" for n in t["core"]]
    L += ["", "## Chapters", "", "| Ch. | Chapter | What it explores |", "|---|---|---|"]
    for c in mine:
        L.append(f"| {c['n']} | [{c['title']}]({rel(src, c['path'])}) | {c['explores']} |")
    L += ["", "Every chapter follows the same arc: **the hook** (a puzzle or a failure), **the big idea**, "
          "**the mathematics**, **the algorithm**, **the experiment** you can run, **the challenge**, and "
          "**loose threads**: the questions the chapter deliberately leaves open."]
    return "\n".join(L) + "\n"


def syllabus(chapters):
    L = ["# Syllabus — Specialization Tracks and Frontier", "",
         "[← Table of Contents](table-of-contents.md) · [Core syllabus](syllabus.md)", "",
         "These chapters extend the 25-chapter core rather than repeat it. Each one links back to the core "
         "chapters it builds on.", ""]
    for t in BOOK["tracks"]:
        L += [f"## Track {t['id']} — {t['title']} (Level {t['level']})", "", t["blurb"], ""]
        for c in [c for c in chapters if c["track"] == t["id"]]:
            line = f"**{c['n']}. {c['title']}** — {c['explores']}"
            if c.get("challenge"):
                line += f" *Challenge: {c['challenge']}*"
            L += [line, ""]
    return "\n".join(L)


def paths(chapters):
    by_n = {c["n"]: c for c in chapters}
    src = "00-front/learning-paths.md"
    L = ["# Learning Paths", "", "[← Table of Contents](table-of-contents.md)", "",
         "You don't need to read the book in order. Each path is a route through the tracks for one goal. "
         "Read the core chapters each track lists first, then follow the path.", ""]
    for p in BOOK["paths"]:
        L += [f"## {p['name']}", "", "Tracks: " + ", ".join(f"Track {x}" for x in p["tracks"]), ""]
        L += [f"{i}. Ch. {n} — [{by_n[n]['title']}]({rel(src, by_n[n]['path'])})" for i, n in enumerate(p["chapters"], 1)]
        L.append("")
    return "\n".join(L)


if __name__ == "__main__":
    core, chapters, apps = core_chapters(), track_chapters(), appendices()
    if "--plan" in sys.argv:
        print(json.dumps(chapters, indent=1, ensure_ascii=False))
        sys.exit()
    write("00-front/table-of-contents.md", toc(core, chapters, apps))
    write("00-front/syllabus-tracks.md", syllabus(chapters))
    write("00-front/learning-paths.md", paths(chapters))
    for t in BOOK["tracks"]:
        write(f"{t['dir']}/README.md", track_index(t, chapters, core))
