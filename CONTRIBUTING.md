# Contributing to The Mathematics Universe

Thanks for helping. A corrected sign, a clearer derivation, a better figure or a whole new chapter all make the book better for every reader.

By taking part you agree to follow the [Code of Conduct](CODE_OF_CONDUCT.md).

## Ways to contribute

| You found or want… | Do this |
|---|---|
| A mathematical error (wrong sign, false statement, broken proof step) | Open an **Erratum** issue, or send a pull request with the fix. Quote the section id, e.g. `§14.7`. |
| A typo, dead link or a graph that doesn't render | Send a pull request straight away; no issue needed. |
| A clearer explanation, extra worked example, better figure | Open a pull request; say briefly why it's clearer. |
| A new chapter, section or track | Open a **Chapter proposal** issue first so we can agree on scope and numbering. |
| A problem with the PDF build, figure scripts or MCP server | Open a **Tooling bug** issue. |

## Setting up

You only need what you're touching.

```bash
git clone https://github.com/The-Studio-of-Sadiq/engineering-mathematics.git
```

- **Reading and editing:** open the folder as an [Obsidian](https://obsidian.md) vault and enable the **Desmos** community plugin. Any Markdown editor works too.
- **Figures:** Python 3.10+ with `numpy`, `matplotlib`, `scipy`.
- **PDF:** Node 18+ and Google Chrome or Microsoft Edge (`cd pdf && npm install`).
- **MCP server:** `pip install -r mcp-server/requirements.txt`.

## House style

Everything is in the [Authoring Guide](00-front/authoring-guide.md). The rules that matter most:

- **Notation** follows [Notation & Conventions](00-front/notation-and-conventions.md): $j$ for the imaginary unit, bold vectors and matrices, cylindrical $(r,\phi,z)$, spherical $(R,\theta,\phi)$.
- **Math:** `$…$` inline and `$$…$$` display, on their own lines. Use `aligned`, never `align`, `\tag`, `\label` or `\newcommand`. Inside tables write `\lvert x\rvert`, not `|x|`.
- **Desmos blocks:** never put `|` inside an equation. Write `\operatorname{abs}(x)`. No double inequalities such as `1<xy<4`.
- **Links:** relative Markdown links (`../part-5-transforms/14-laplace-transforms.md#147-…`), not `[[wikilinks]]`.
- **Core chapters (1–25)** are thorough reference chapters. **Track chapters (26–99)** are shorter and follow the seven-part arc: hook, big idea, mathematics, algorithm, experiment, challenge, loose threads. They plant `[!mystery]` clues that are answered elsewhere, and their challenges have **no answer key**. Please don't add solutions.
- **Figures** are drawn by a script `figures/src/chNN.py` using `figures/style.py`, never pasted in from elsewhere. Commit the script together with the SVG it writes.
- **Check every number.** If a worked example says 0.8433, a Python one-liner should agree.

## Before you open a pull request

Run the checks; CI runs the same ones.

```bash
python tools/check_book.py           # math delimiters, links, headings, stray characters
python tools/fix_desmos.py --check   # Desmos syntax the Obsidian plugin can't parse
python tools/gen_nav.py              # regenerate TOC and track pages (then commit any changes)
```

If you changed a figure script:

```bash
python figures/build_all.py chNN
```

If you changed anything that affects the PDF and want to check it, `cd pdf && node build.mjs --only tracks` (or `core`, `appendix`) builds a partial PDF. You don't need to commit `The-Mathematics-Universe.pdf`; maintainers rebuild it for releases.

## Adding a chapter

1. Agree on the number and slug in a Chapter proposal issue.
2. Add the chapter to [`tools/book.json`](tools/book.json) (track chapters) and run `python tools/gen_nav.py`.
3. Write the chapter following the guide, plus `figures/src/chNN.py`.
4. Run the checks above.

## Commits and pull requests

- Keep each pull request to one topic (one chapter, one fix, one figure).
- Write commit messages in the imperative: "Fix sign in §14.7 final-value example".
- Fill in the pull-request checklist.

## Licensing of contributions

By contributing you agree that your contribution is licensed under the repository's licenses: **CC BY 4.0** for book content ([LICENSE-CONTENT.md](LICENSE-CONTENT.md)) and **MIT** for code ([LICENSE](LICENSE)). Only contribute material you wrote yourself or that is compatibly licensed, and credit any source you adapt.
