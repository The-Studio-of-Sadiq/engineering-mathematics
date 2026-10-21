# Changelog

All notable changes to *The Mathematics Universe*. Dates are ISO 8601.

## [1.0.0] — 2026-10-22

First complete edition.

### Added
- **Level 1 — Core:** 25 chapters of engineering mathematics (Parts I–VII), 13,000–21,000 words each, with 245 Desmos graphs.
- **Appendix A–E:** complete proofs of 61 central results.
- **Level 2 — Specializations:** Tracks A–I, chapters 26–88 and 99 (algorithms; probability and information; game theory; optimization; abstract mathematics; logic and foundations; graphics and simulation; cryptography and coding; dynamical systems and complexity). Each chapter has a hook, big idea, mathematics, algorithm, experiment, challenge and loose threads, with `[!mystery]` missing clues.
- **Level 3 — Research Frontier:** Track J, chapters 89–98.
- 198 original figures drawn by scripts in `figures/`, track roadmaps and a universe map.
- `The-Mathematics-Universe.pdf`: 1,422 pages with bookmarks, cross-reference links, KaTeX math and Desmos-rendered graphs, built by `pdf/build.mjs`.
- MCP server (`mcp-server/`) for searching and reading the book from Claude and other MCP clients.
- Tooling: `tools/book.json` (structure), `tools/gen_nav.py` (navigation), `tools/fix_desmos.py` (Desmos syntax), `tools/check_book.py` (consistency checks), CI workflow.
- Community files: licenses (CC BY 4.0 for content, MIT for code), contributing guide, code of conduct, security policy, citation file, issue and pull-request templates.

### Fixed
- 29 Desmos graphs in the core and appendix that the Obsidian plugin could not parse (absolute values written with `|`, double inequalities, implicit list products).
