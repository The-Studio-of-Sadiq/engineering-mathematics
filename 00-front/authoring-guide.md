# Authoring Guide

[← Table of Contents](table-of-contents.md)

Rules every chapter follows, so the book reads as one piece.

## File layout

- One chapter = one file: `part-N-<name>/NN-<slug>.md` (slugs fixed in the [Table of Contents](table-of-contents.md)).
- Plain Markdown + LaTeX. No images to maintain; graphs are `desmos-graph` blocks.
- Links are relative Markdown links (`[Ch. 14](../part-5-transforms/14-laplace-transforms.md)`), not `[[wikilinks]]`, so they work on GitHub and in Obsidian.

Appendix files are `appendix/X-<slug>.md` (X = A, B, …), with sections `## X.k Result name`. Each holds a theorem statement, a full proof ending in ∎, and remarks.

## Chapter skeleton

```markdown
# Chapter N — Title

[← Previous](…) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](…)

> **Prerequisites:** Ch. X, Ch. Y
> **Learning outcomes:** After this chapter you can …  (3–6 bullets)

## N.1 Section title          ← one section per syllabus bullet, in syllabus order
…definitions, key results, a worked example where it helps…

## N.k Section title
…

## Summary                    ← key formulas table
## Further reading            ← 2–4 standard texts
```

No exercises, problems, or answers sections — this is a reference text; worked examples inside sections carry the teaching.

## Math

- Inline `$…$`, display `$$…$$` on their own lines with a blank line before and after.
- Imaginary unit is $j$. Use `\mathbf{}` for vectors/matrices, `\operatorname{}` for named operators.
- Avoid `\begin{align}` numbering macros; use `\begin{aligned}…\end{aligned}` inside `$$…$$` (works in both KaTeX and MathJax).
- Do not use `\newcommand`, `\label`, `\eqref`, `\tag` — they are not portable.
- Inside tables, escape the pipe: write `\lvert x\rvert` instead of `|x|`.

## Callouts

Use `> [!definition]`, `> [!theorem]`, `> [!example]`, `> [!tip]` (engineering note), `> [!warning]` (pitfall). Every line of a callout starts with `>`. Display math inside a callout:

```markdown
> [!theorem] Euler's formula
> $$
> e^{j\theta} = \cos\theta + j\sin\theta
> $$
```

## Graphs (`desmos-graph`)

Rendered by the [obsidian-desmos](https://github.com/Nigecat/obsidian-desmos) plugin. Format: optional settings, a `---` line, then one equation per line.

````markdown
```desmos-graph
left=-1; right=10; bottom=-1.2; top=1.2
---
y=e^{-0.5x}\cos(3x)|blue
y=e^{-0.5x}|red|dashed
y=-e^{-0.5x}|red|dashed
(2,0.5)|label:point
```
````

Rules that keep graphs rendering:

- Settings: only `left`, `right`, `top`, `bottom`, `width`, `height`, `grid`, `degreeMode`. Separate with `;`.
- Each equation is Desmos LaTeX: `\sin(x)`, `\frac{a}{b}`, `e^{x}`, `\sqrt{x}`, `\left|x\right|`, `\sum_{n=1}^{9}`.
- Style options after `|`: a colour (`red`, `green`, `blue`, `yellow`, `orange`, `purple`, `grey`, `black` or `#rrggbb`), `dashed`/`dotted`, a domain restriction such as `0<x<5`, `label:text`.
- Variables are `x` and `y` (or `r`, `\theta` for polar). Fix constants numerically — do not rely on sliders.
- Use graphs where seeing the shape teaches something (2–5 per chapter). Put a one-line caption in italics under each graph.
- Desmos cannot draw 3-D surfaces, vector fields or histograms well; describe those in words or tables instead.

## Depth

Each chapter is a full textbook chapter: **30–40 printed pages, i.e. roughly 12,000–16,000 words (about 1,800–2,800 Markdown lines)**. Every syllabus section gets:

- motivation (why an engineer needs it),
- precise definitions and theorem statements,
- derivations or proof sketches for the central results (enough that the reader could reconstruct them),
- 2–4 worked examples of increasing difficulty, each solved step by step and checked numerically,
- a pitfall and/or engineering note,
- graphs where shape matters.

## Tone and scope

- Textbook-standard results, explained thoroughly; derive what is central, cite what is routine.
- EE examples (RC/RLC circuits, phasors, signals, fields, noise) motivate methods but the chapter stays a math chapter (see book boundary in the [README](../README.md)).
- Check every equation for dimensional consistency and every worked example numerically.
