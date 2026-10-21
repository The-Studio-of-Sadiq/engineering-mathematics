"""MCP server exposing The Mathematics Universe book (core, tracks and appendix).

Read-only tools for listing chapters, reading sections, searching the text,
and pulling formula summaries, notation and Desmos graphs. Chapters are
addressed by number ("14"), appendices by letter ("D"); sections by their
heading number ("14.6", "D.3").

Run: python mcp-server/server.py   (stdio transport)
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError

BOOK_ROOT = Path(__file__).resolve().parent.parent

HEADING_RE = re.compile(r"^(#{1,6})\s+(.*\S)\s*$")
SECTION_ID_RE = re.compile(r"^([0-9]+|[A-Z])\.(\d+)\b")
DESMOS_RE = re.compile(r"```desmos-graph\n(.*?)```[ \t]*\n(?:\s*\n)*(\*[^\n]+\*)?", re.S)


@dataclass
class Section:
    id: str | None  # "14.6", or None for unnumbered (Summary, Further reading, preamble)
    title: str
    start: int  # line index of the heading
    end: int  # exclusive


@dataclass
class Unit:
    key: str  # "14" or "D"
    kind: str  # "chapter" or "appendix"
    title: str
    part: str
    path: Path
    lines: list[str] = field(default_factory=list)
    sections: list[Section] = field(default_factory=list)

    @property
    def rel(self) -> str:
        return self.path.relative_to(BOOK_ROOT).as_posix()


def _parse(unit: Unit) -> Unit:
    unit.lines = unit.path.read_text(encoding="utf-8").splitlines()
    in_code = False
    heads: list[tuple[int, str]] = []
    for i, line in enumerate(unit.lines):
        if line.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        m = HEADING_RE.match(line)
        if not m:
            continue
        level, text = len(m.group(1)), m.group(2)
        if level == 1 and not unit.title:
            unit.title = text
        elif level == 2:
            heads.append((i, text))
    bounds = [i for i, _ in heads] + [len(unit.lines)]
    for (start, text), end in zip(heads, bounds[1:]):
        sid = SECTION_ID_RE.match(text)
        unit.sections.append(Section(sid.group(0) if sid else None, text, start, end))
    return unit


def _load() -> dict[str, Unit]:
    units: dict[str, Unit] = {}
    for path in sorted(BOOK_ROOT.glob("part-*/*.md")):
        m = re.match(r"(\d+)-", path.name)
        if m:
            key = str(int(m.group(1)))
            units[key] = _parse(Unit(key, "chapter", "", path.parent.name, path))
    for path in sorted(BOOK_ROOT.glob("tracks/*/*.md")):
        m = re.match(r"(\d+)-", path.name)
        if m:
            key = str(int(m.group(1)))
            units[key] = _parse(Unit(key, "track chapter", "", path.parent.name, path))
    for path in sorted((BOOK_ROOT / "appendix").glob("*.md")):
        m = re.match(r"([A-Z])-", path.name)
        if m:
            units[m.group(1)] = _parse(Unit(m.group(1), "appendix", "", "appendix", path))
    return units


def _units() -> dict[str, Unit]:
    # Reload on every call so edits in Obsidian are visible without restarting.
    return _load()


def _get(chapter: str) -> Unit:
    key = chapter.strip().upper().removeprefix("CH").removeprefix("APPENDIX").strip(" .")
    key = str(int(key)) if key.isdigit() else key
    units = _units()
    if key not in units:
        raise ToolError(f"Unknown chapter/appendix {chapter!r}. Known: {', '.join(units)}")
    return units[key]


def _text(unit: Unit, sec: Section) -> str:
    return "\n".join(unit.lines[sec.start : sec.end]).strip()


def _section_for_line(unit: Unit, line_no: int) -> Section | None:
    for sec in unit.sections:
        if sec.start <= line_no < sec.end:
            return sec
    return None


mcp = MCPServer(
    name="engineering-math-book",
    title="The Mathematics Universe",
    instructions=(
        "The Mathematics Universe: a 25-chapter engineering-mathematics core (keys '1'..'25'), "
        "specialization and frontier track chapters (keys '26'..'99': algorithms, probability and "
        "information, game theory, optimization, abstract math, logic, graphics, cryptography, "
        "dynamical systems, frontier topics) and appendices of proofs (keys 'A'..'E'). "
        "Use list_chapters or search first, then get_section to read exact text. "
        "Math is LaTeX ($...$, $$...$$); imaginary unit is j."
    ),
)


@mcp.tool()
def list_chapters() -> list[dict]:
    """List every chapter and appendix with its key, title, part, file path and section count."""
    return [
        {"key": u.key, "kind": u.kind, "title": u.title, "part": u.part, "path": u.rel,
         "sections": sum(1 for s in u.sections if s.id), "words": sum(len(l.split()) for l in u.lines)}
        for u in _units().values()
    ]


@mcp.tool()
def get_outline(chapter: str) -> dict:
    """Return the section headings of a chapter ('14') or appendix ('D'), including ### subheadings."""
    unit = _get(chapter)
    outline = []
    for sec in unit.sections:
        subs, in_code = [], False
        for line in unit.lines[sec.start + 1 : sec.end]:
            if line.startswith("```"):
                in_code = not in_code
            elif not in_code and line.startswith("### "):
                subs.append(line[4:].strip())
        outline.append({"id": sec.id, "title": sec.title, "subheadings": subs})
    return {"key": unit.key, "title": unit.title, "path": unit.rel, "sections": outline}


@mcp.tool()
def get_section(section: str, max_chars: int = 40000) -> str:
    """Return the full Markdown of one section, e.g. '14.6' or 'D.3'.

    Also accepts '<chapter>:<title words>' (e.g. '14:Summary', '7:Further reading') for
    unnumbered sections. Output is truncated to max_chars.
    """
    if ":" in section:
        chapter, name = section.split(":", 1)
        unit = _get(chapter)
        matches = [s for s in unit.sections if name.strip().lower() in s.title.lower()]
    else:
        chapter = section.split(".", 1)[0]
        unit = _get(chapter)
        want = f"{unit.key}.{section.split('.', 1)[1]}" if "." in section else None
        matches = [s for s in unit.sections if s.id == want]
    if not matches:
        ids = ", ".join(s.id or s.title for s in unit.sections)
        raise ToolError(f"No section {section!r} in {unit.title}. Available: {ids}")
    text = _text(unit, matches[0])
    if len(text) > max_chars:
        text = text[:max_chars] + f"\n\n[... truncated at {max_chars} chars of {len(text)}]"
    return f"<!-- {unit.rel} -->\n{text}"


@mcp.tool()
def get_summary(chapter: str) -> str:
    """Return a chapter's Summary section (its table of key formulas)."""
    return get_section(f"{chapter}:Summary")


@mcp.tool()
def search(query: str, chapter: str | None = None, limit: int = 15) -> list[dict]:
    """Case-insensitive search over the book. All words in the query must appear on a line
    (or a quoted "exact phrase" must match); hyphens count as spaces. Returns one hit per
    section (its best-matching line), ranked with section-title and heading matches first,
    then by how often the query words occur in the section."""
    norm = lambda s: s.lower().replace("-", " ")
    phrase = re.fullmatch(r'"(.+)"', query.strip())
    terms = [norm(phrase.group(1))] if phrase else norm(query).split()
    if not terms:
        return []
    # Match at word starts so "rms" does not hit "transforms" but "theorem" hits "theorems".
    pats = [re.compile(r"(?<![a-z0-9])" + re.escape(t)) for t in terms]
    has_all = lambda s: all(p.search(s) for p in pats)
    count = lambda s: sum(len(p.findall(s)) for p in pats)
    units = [_get(chapter)] if chapter else list(_units().values())
    best: dict[tuple[str, int], dict] = {}
    for unit in units:
        for i, line in enumerate(unit.lines):
            low = norm(line)
            if not has_all(low):
                continue
            sec = _section_for_line(unit, i)
            sec_text = norm(_text(unit, sec)) if sec else low
            title_hit = sec is not None and has_all(norm(sec.title))
            density = count(sec_text) / max(1, len(sec_text) / 5000)
            score = 100 * title_hit + 20 * line.startswith("#") + 10 * ("[!theorem]" in line) + density
            snippet = line.strip()
            if len(snippet) > 300:
                pos = pats[0].search(low).start()
                snippet = "…" + snippet[max(0, pos - 120) : pos + 180] + "…"
            key = (unit.key, sec.start if sec else -1)
            if key not in best or score > best[key]["score"]:
                best[key] = {"chapter": unit.key, "chapter_title": unit.title,
                             "section": sec.id if sec else None, "section_title": sec.title if sec else None,
                             "line": i + 1, "snippet": snippet, "score": round(score, 1)}
    return sorted(best.values(), key=lambda h: -h["score"])[:limit]


@mcp.tool()
def get_graphs(chapter: str) -> list[dict]:
    """Return the Desmos graphs of a chapter/appendix: section, settings, equations and caption.
    Equations are Desmos LaTeX and can be pasted line by line into desmos.com/calculator."""
    unit = _get(chapter)
    text = "\n".join(unit.lines)
    graphs = []
    for m in DESMOS_RE.finditer(text):
        body = m.group(1)
        settings, _, eqs = body.partition("\n---\n") if "\n---\n" in body else ("", "", body)
        line_no = text.count("\n", 0, m.start())
        sec = _section_for_line(unit, line_no)
        graphs.append({"section": sec.id if sec else None, "line": line_no + 1,
                       "settings": settings.strip(),
                       "equations": [e for e in eqs.strip().splitlines() if e.strip()],
                       "caption": (m.group(2) or "").strip("* ")})
    return graphs


@mcp.tool()
def get_notation() -> str:
    """Return the book's notation and conventions page (symbols, coordinate systems, units)."""
    return (BOOK_ROOT / "00-front" / "notation-and-conventions.md").read_text(encoding="utf-8")


@mcp.resource("book://toc", mime_type="text/markdown", description="Table of contents")
def toc() -> str:
    return (BOOK_ROOT / "00-front" / "table-of-contents.md").read_text(encoding="utf-8")


@mcp.prompt()
def explain_topic(topic: str) -> str:
    """Explain a topic using the book as the source of truth."""
    return (
        f"Explain '{topic}' for an electrical engineering student. First call search with the key "
        "terms, then get_section on the most relevant sections (and the matching appendix proof if "
        "one exists). Base the explanation on the book's text and notation, cite section ids, and "
        "include one worked example from the book."
    )


if __name__ == "__main__":
    mcp.run()
