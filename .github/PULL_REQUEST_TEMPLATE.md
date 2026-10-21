## What this changes

<!-- One topic per pull request: a fix, a figure, a section, a chapter. Name the section ids, e.g. §14.7. -->

## Why

<!-- For a correction: the derivation, counterexample or numerical check that shows the old text was wrong. -->

## Checklist

- [ ] `python tools/check_book.py` passes
- [ ] `python tools/fix_desmos.py --check` reports 0
- [ ] Ran `python tools/gen_nav.py` if chapters or `tools/book.json` changed, and committed the result
- [ ] New or changed figures come from a script in `figures/src/`, committed with its SVG
- [ ] Every number I added was checked (e.g. with Python)
- [ ] Notation and style follow the [Authoring Guide](../00-front/authoring-guide.md)
- [ ] No solutions to track-chapter challenges or `[!mystery]` questions were added
- [ ] I agree my contribution is licensed under CC BY 4.0 (content) / MIT (code)
