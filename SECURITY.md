# Security Policy

This repository is mostly a book. The parts that run code are:

- the **MCP server** (`mcp-server/`), which reads the Markdown files on the local machine and serves them to MCP clients;
- the **PDF build** (`pdf/`), which drives a local headless Chrome and loads the Desmos API;
- the **figure scripts** (`figures/`) and **tools** (`tools/`).

## Supported versions

Only the latest commit on `main` is supported.

## Reporting a vulnerability

Please **do not** open a public issue for security problems. Email **princekhan7748@gmail.com** with:

- what the problem is and which component it affects;
- steps to reproduce, or a proof of concept;
- the impact you expect.

You should get a reply within a week. Once a fix is ready we will credit you, if you wish, in the release notes.

## A note on the cryptography chapters

Track H (chapters 77–82) teaches the mathematics of cryptography with **deliberately tiny, insecure toy parameters**. Those code snippets are teaching examples, not security tools. A weakness in a toy example is not a vulnerability, but a mathematical error in the text is welcome as an Erratum issue.
