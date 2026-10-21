# Book MCP Server

An [MCP](https://modelcontextprotocol.io) server that lets Claude (or any MCP client) read this book: list chapters, pull exact sections, search the text, and fetch formula summaries, notation and Desmos graphs. It is read-only and re-reads the Markdown on every call, so edits made in Obsidian show up immediately.

## Setup

```bash
pip install -r mcp-server/requirements.txt
```

**Claude Code:** the repo's [`.mcp.json`](../.mcp.json) registers the server automatically. Open Claude Code in the repo root and approve `engineering-math-book` when prompted. It runs `python mcp-server/server.py` from the repo root, so `python` must be the interpreter that has `mcp` installed.

**Claude Desktop / other clients:** add to the client's MCP config, using absolute paths:

```json
{
  "mcpServers": {
    "engineering-math-book": {
      "command": "python",
      "args": ["C:/obsidian/EM/mcp-server/server.py"]
    }
  }
}
```

## Addressing

- Chapters by number: `"14"`. Appendices by letter: `"D"`.
- Numbered sections by heading id: `"14.7"`, `"D.3"`.
- Unnumbered sections by `<chapter>:<title words>`: `"14:Summary"`, `"7:Further reading"`.

## Tools

| Tool | Arguments | Returns |
|---|---|---|
| `list_chapters` | — | key, kind, title, part, path, section and word counts for every chapter and appendix |
| `get_outline` | `chapter` | `##` sections with their `###` subheadings |
| `get_section` | `section`, `max_chars=40000` | full Markdown of one section |
| `get_summary` | `chapter` | the chapter's key-formula table |
| `search` | `query`, `chapter=None`, `limit=15` | one hit per section: id, title, line, snippet, score. All words must match at word starts; hyphens count as spaces; `"quoted phrase"` for exact phrases. Section-title and heading matches rank first. |
| `get_graphs` | `chapter` | each Desmos block's section, settings, equations and caption |
| `get_notation` | — | the notation and conventions page |

Also: resource `book://toc` (table of contents) and prompt `explain_topic(topic)`, which tells the model to answer from the book and cite section ids.

## How it parses

A chapter is any `part-*/NN-*.md` or `tracks/*/NN-*.md`; an appendix is any `appendix/X-*.md`. The `#` heading is the title. `##` headings split sections, and a leading `14.7` or `D.3` in the heading becomes the section id. Headings inside code fences are ignored. Follow the [authoring guide](../00-front/authoring-guide.md) and new content is picked up without code changes.
