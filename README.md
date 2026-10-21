# Book 1 — Engineering Mathematics Core

> **Objective.** Develop the mathematical foundation needed to understand, model, analyze, and simulate electrical and electronic systems.

**Boundary.** This book teaches *mathematical methods*. Applications such as Fourier-based filtering, circuit simulation, and motor modeling are treated in their relevant core books or in *Applications & Integrations*. Engineering examples appear here only to motivate and illustrate a method.

---

## How to read this repo

- Every chapter is a single Markdown file with LaTeX math (`$…$` inline, `$$…$$` display). It renders in **Obsidian**, **GitHub**, **VS Code** (Markdown preview), and any KaTeX/MathJax renderer.
- Graphs are written as `desmos-graph` code blocks — plain text equations, no image files to maintain. See [Graphs](#graphs) below.
- Start with [Notation & Conventions](00-front/notation-and-conventions.md), then follow the [Table of Contents](00-front/table-of-contents.md).
- Contributors: read the [Authoring Guide](00-front/authoring-guide.md) before adding content.

## Graphs

Graphs use the [Obsidian Desmos plugin](https://github.com/Nigecat/obsidian-desmos), which renders the Desmos graphing calculator from a fenced block:

````markdown
```desmos-graph
left=-7; right=7; top=2; bottom=-2
---
y=\sin(x)|blue
y=\cos(x)|red|dashed
```
````

**Setup (Obsidian):** Settings → Community plugins → Browse → search **"Desmos"** → Install → Enable.

Outside Obsidian (e.g. GitHub), the block shows as readable source — the equations themselves — which you can paste into [desmos.com/calculator](https://www.desmos.com/calculator) line by line.

## Structure

| Part | Folder | Chapters |
|---|---|---|
| I — Mathematical Foundations | [part-1-foundations](part-1-foundations/) | 1–3 |
| II — Calculus | [part-2-calculus](part-2-calculus/) | 4–7 |
| III — Linear Algebra | [part-3-linear-algebra](part-3-linear-algebra/) | 8–10 |
| IV — Differential Equations | [part-4-differential-equations](part-4-differential-equations/) | 11–13 |
| V — Transforms and Mathematical Analysis | [part-5-transforms](part-5-transforms/) | 14–18 |
| VI — Probability, Statistics, and Stochastic Mathematics | [part-6-probability-statistics](part-6-probability-statistics/) | 19–21 |
| VII — Numerical and Computational Mathematics | [part-7-numerical-computational](part-7-numerical-computational/) | 22–25 |
| Appendix — Derivations and Proofs | [appendix](appendix/) | A–E |

Full chapter list: [Table of Contents](00-front/table-of-contents.md).

The chapters state key results with short justifications. The appendix gives complete proofs, with precise hypotheses and counterexamples, and links back to where each result is used.

## MCP server

[`mcp-server/`](mcp-server/) lets Claude and other MCP clients read the book: list chapters, fetch exact sections, search the text, and pull formula summaries and graphs. Claude Code picks it up from [`.mcp.json`](.mcp.json). For setup and the tool list, see the [server README](mcp-server/README.md).
