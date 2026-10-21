// Build The Mathematics Universe as one PDF.
//
//   cd pdf && npm install && node build.mjs            -> ../The-Mathematics-Universe.pdf
//   node build.mjs --desmos-only                       -> only (re)render Desmos graphs into cache/
//   node build.mjs --only core|tracks|appendix         -> a partial PDF (faster, for checking)
//
// Pipeline: collect chapters in reading order (tools/book.json + core folders), replace each
// ```desmos-graph block with a PNG rendered by the real Desmos calculator in headless Chrome,
// turn Obsidian callouts into styled boxes, render math with KaTeX, rewrite .md links into
// in-document anchors, and print the whole book with Chrome (with a PDF outline).
//
// Needs Chrome or Edge installed. Desmos needs network access and an API key: set DESMOS_API_KEY,
// otherwise the public demo key from the Desmos API docs is used.

import { createHash } from "node:crypto";
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import katex from "katex";
import MarkdownIt from "markdown-it";
import { PDFDocument } from "pdf-lib";
import { chromium } from "playwright-core";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");
const CACHE = path.join(HERE, "cache", "desmos");
const OUT_HTML = path.join(HERE, "cache", "book.html");
const BOOK = JSON.parse(fs.readFileSync(path.join(ROOT, "tools", "book.json"), "utf8"));
const args = process.argv.slice(2);
const ONLY = args.includes("--only") ? args[args.indexOf("--only") + 1] : null;
const OUT_PDF = path.join(ROOT, ONLY ? `The-Mathematics-Universe-${ONLY}.pdf` : "The-Mathematics-Universe.pdf");
const DESMOS_KEY = process.env.DESMOS_API_KEY || "dcb31709b452b1cf9dc26972add0fda6";

// ---------------------------------------------------------------- reading order

const CORE_PARTS = [
  ["Part I — Mathematical Foundations", "part-1-foundations"],
  ["Part II — Calculus", "part-2-calculus"],
  ["Part III — Linear Algebra", "part-3-linear-algebra"],
  ["Part IV — Differential Equations", "part-4-differential-equations"],
  ["Part V — Transforms and Mathematical Analysis", "part-5-transforms"],
  ["Part VI — Probability, Statistics, and Stochastic Mathematics", "part-6-probability-statistics"],
  ["Part VII — Numerical and Computational Mathematics", "part-7-numerical-computational"],
];

const mdIn = (dir, re) =>
  fs.existsSync(path.join(ROOT, dir))
    ? fs.readdirSync(path.join(ROOT, dir)).filter((f) => re.test(f)).sort().map((f) => `${dir}/${f}`)
    : [];

/** [{divider?: string, file?: string}] in reading order. */
function readingOrder() {
  const items = [{ file: "00-front/notation-and-conventions.md" }];
  if (!ONLY || ONLY === "core") {
    items.push({ divider: "Level 1 — Core", blurb: BOOK.levels[0].blurb });
    for (const [part, dir] of CORE_PARTS) {
      items.push({ divider: part, minor: true });
      for (const f of mdIn(dir, /^\d\d-.*\.md$/)) items.push({ file: f });
    }
  }
  if (!ONLY || ONLY === "tracks") {
    for (const level of [2, 3]) {
      const L = BOOK.levels.find((l) => l.level === level);
      items.push({ divider: `Level ${level} — ${L.name}`, blurb: L.blurb });
      for (const t of BOOK.tracks.filter((t) => t.level === level)) {
        if (fs.existsSync(path.join(ROOT, t.dir, "README.md"))) items.push({ file: `${t.dir}/README.md` });
        for (const c of t.chapters) {
          const f = `${t.dir}/${String(c.n).padStart(2, "0")}-${c.slug}.md`;
          if (fs.existsSync(path.join(ROOT, f))) items.push({ file: f });
          else console.warn(`missing chapter ${f}`);
        }
      }
    }
  }
  if (!ONLY || ONLY === "appendix") {
    items.push({ divider: "Appendix — Derivations and Proofs" });
    for (const f of mdIn("appendix", /^[A-Z]-.*\.md$/)) items.push({ file: f });
  }
  return items;
}

// ---------------------------------------------------------------- helpers

const docId = (rel) => "d-" + rel.replace(/\.md$/, "").replace(/[^A-Za-z0-9]+/g, "-").toLowerCase();

/** GitHub/Obsidian-style heading slug, so existing #fragments in links resolve. */
const slugify = (s) =>
  s.toLowerCase().trim()
    .replace(/<[^>]+>/g, "")
    .replace(/[^\p{L}\p{N}\s_-]/gu, "")
    .replace(/\s/g, "-");

const escapeHtml = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");

// ---------------------------------------------------------------- Desmos

const COLORS = {
  red: "#c74440", blue: "#2d70b3", green: "#388c46", purple: "#6042a6", orange: "#fa7e19",
  black: "#000000", grey: "#6b6b6b", gray: "#6b6b6b", yellow: "#d4a017", magenta: "#c837c4",
  cyan: "#00a3c4", white: "#ffffff",
};

/** Parse an obsidian-desmos block into {settings, expressions}. */
function parseDesmos(body) {
  let settingsText = "", eqText = body;
  const sep = body.split(/\n-{3,}\s*\n/);
  if (sep.length > 1) [settingsText, eqText] = [sep[0], sep.slice(1).join("\n")];
  const settings = {};
  for (const kv of settingsText.split(/[;\n]/)) {
    const m = kv.match(/^\s*([A-Za-z]+)\s*=\s*(.+?)\s*$/);
    if (m) settings[m[1]] = m[2];
  }
  const expressions = [];
  for (const raw of eqText.split("\n")) {
    if (!raw.trim()) continue;
    const [latex, ...opts] = raw.split("|");
    const e = { latex: latex.trim(), restrictions: [] };
    for (let o of opts) {
      o = o.trim();
      if (!o) continue;
      const lo = o.toLowerCase();
      if (COLORS[lo]) e.color = COLORS[lo];
      else if (/^#[0-9a-f]{3,8}$/i.test(o)) e.color = o;
      else if (["dashed", "dotted", "solid"].includes(lo)) e.lineStyle = lo.toUpperCase();
      else if (lo === "hidden") e.hidden = true;
      else if (lo.startsWith("label:")) e.label = o.slice(6).trim();
      else if (lo === "label") e.label = "";
      else if (/^[a-z]+$/.test(lo)) { /* unknown keyword (open, cross, noline, ...) */ }
      else e.restrictions.push(o);
    }
    expressions.push(e);
  }
  return { settings, expressions };
}

const desmosHash = (body) => createHash("sha1").update(body.trim()).digest("hex").slice(0, 16);
const blockFiles = new Map(); // hash -> Set of source files

/** Render every not-yet-cached Desmos block to cache/desmos/<hash>.png. */
async function renderDesmos(blocks, browser) {
  fs.mkdirSync(CACHE, { recursive: true });
  const todo = [...blocks.entries()].filter(([h]) => !fs.existsSync(path.join(CACHE, `${h}.png`)));
  console.log(`desmos: ${blocks.size} graphs, ${todo.length} to render`);
  if (!todo.length) return;
  const page = await browser.newPage({ viewport: { width: 900, height: 600 } });
  await page.setContent(`<!doctype html><html><body style="margin:0">
    <div id="calc" style="width:720px;height:432px"></div>
    <script src="https://www.desmos.com/api/v1.10/calculator.js?apiKey=${DESMOS_KEY}"></script>
    </body></html>`, { waitUntil: "load", timeout: 120000 });
  await page.waitForFunction(() => window.Desmos, null, { timeout: 120000 });
  await page.evaluate(() => {
    window.calc = Desmos.GraphingCalculator(document.getElementById("calc"), {
      expressions: false, settingsMenu: false, zoomButtons: false, keypad: false, border: false,
      lockViewport: true, showResetButtonOnGraphpaper: false,
    });
  });
  let n = 0, failed = 0;
  for (const [hash, body] of todo) {
    const spec = parseDesmos(body);
    const result = await page.evaluate(async (spec) => {
      const s = spec.settings, calc = window.calc;
      calc.setBlank();
      const num = (v, d) => (v !== undefined && !isNaN(parseFloat(v)) ? parseFloat(v) : d);
      calc.updateSettings({
        degreeMode: /deg/i.test(s.degreeMode || ""),
        showGrid: !/false/i.test(s.grid || ""),
        showXAxis: true, showYAxis: true,
        xAxisNumbers: !/true/i.test(s.hideAxisNumbers || ""), yAxisNumbers: !/true/i.test(s.hideAxisNumbers || ""),
        xAxisLabel: s.xAxisLabel || "", yAxisLabel: s.yAxisLabel || "",
        xAxisScale: /true/i.test(s.xAxisLogarithmic || "") ? "logarithmic" : "linear",
        yAxisScale: /true/i.test(s.yAxisLogarithmic || "") ? "logarithmic" : "linear",
      });
      calc.setMathBounds({ left: num(s.left, -10), right: num(s.right, 10), bottom: num(s.bottom, -7), top: num(s.top, 7) });
      spec.expressions.forEach((e, i) => {
        let latex = e.latex;
        for (const r of e.restrictions) latex += `\\left\\{${r}\\right\\}`;
        const expr = { id: "e" + i, latex };
        if (e.color) expr.color = e.color;
        if (e.lineStyle) expr.lineStyle = Desmos.Styles[e.lineStyle];
        if (e.hidden) expr.hidden = true;
        if (e.label !== undefined) { expr.showLabel = true; if (e.label) expr.label = e.label; }
        calc.setExpression(expr);
      });
      await new Promise((r) => setTimeout(r, 60));
      const errors = spec.expressions
        .map((e, i) => ({ latex: e.latex, a: calc.expressionAnalysis["e" + i] }))
        .filter((x) => x.a && x.a.isError)
        .map((x) => ({ latex: x.latex, error: x.a.errorMessage }));
      const png = await new Promise((resolve) =>
        calc.asyncScreenshot({ width: 720, height: 432, targetPixelRatio: 1.5, mode: "contain" }, resolve));
      return { png, errors };
    }, spec);
    fs.writeFileSync(path.join(CACHE, `${hash}.png`), Buffer.from(result.png.split(",")[1], "base64"));
    const errFile = path.join(CACHE, `${hash}.err.json`);
    if (result.errors.length) {
      fs.writeFileSync(errFile, JSON.stringify({ files: [...(blockFiles.get(hash) || [])], errors: result.errors }, null, 1));
      failed++;
    } else if (fs.existsSync(errFile)) fs.rmSync(errFile);
    if (++n % 50 === 0) console.log(`  rendered ${n}/${todo.length}`);
  }
  console.log(`desmos: rendered ${n}, ${failed} with expression errors`);
  await page.close();
}

/** Print every cached Desmos error that belongs to a graph still in the book. */
function reportDesmosErrors(blocks) {
  const lines = [];
  for (const h of blocks.keys()) {
    const f = path.join(CACHE, `${h}.err.json`);
    if (!fs.existsSync(f)) continue;
    const { errors } = JSON.parse(fs.readFileSync(f, "utf8"));
    for (const e of errors) lines.push(`${[...(blockFiles.get(h) || [])].join(", ")}: ${e.latex}  -->  ${e.error}`);
  }
  fs.writeFileSync(path.join(HERE, "cache", "desmos-errors.txt"), lines.join("\n") + (lines.length ? "\n" : ""));
  if (lines.length) console.log(`desmos: ${lines.length} expression errors, see pdf/cache/desmos-errors.txt`);
}

// ---------------------------------------------------------------- Markdown -> HTML

const md = new MarkdownIt({ html: true, linkify: false, typographer: false });

/** Protect math from markdown-it: swap it for placeholder tokens, rendered by KaTeX afterwards. */
function extractMath(text, store) {
  const put = (tex, display) => {
    store.push({ tex, display });
    return `KTXPH${store.length - 1}ENDPH`;
  };
  // Split around fenced code and inline code so their dollars are left alone.
  const parts = text.split(/(```[\s\S]*?```|`[^`\n]+`)/g);
  return parts.map((p, i) => {
    if (i % 2) return p;
    p = p.replace(/\$\$([\s\S]+?)\$\$/g, (_, tex) => `\n\n${put(tex, true)}\n\n`);
    p = p.replace(/(^|[^\\$])\$(?!\s)([^$\n]+?)(?<![\s\\])\$(?!\d)/g, (_, pre, tex) => pre + put(tex, false));
    return p;
  }).join("");
}

function renderMath(html, store, errors, file) {
  return html.replace(/KTXPH(\d+)ENDPH/g, (_, i) => {
    const { tex, display } = store[+i];
    try {
      return katex.renderToString(tex.trim(), { displayMode: display, throwOnError: true, strict: "ignore", output: "html", trust: false });
    } catch (e) {
      errors.push(`${file}: ${e.message.split("\n")[0]}`);
      return katex.renderToString(tex.trim(), { displayMode: display, throwOnError: false, strict: "ignore", output: "html" });
    }
  });
}

/** Obsidian callouts (> [!type] Title ...) -> <div class="callout type"> with the body as Markdown. */
function convertCallouts(text) {
  const lines = text.split("\n"), out = [];
  let inFence = false;
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (/^\s*```/.test(line)) inFence = !inFence;
    const m = !inFence && line.match(/^>\s*\[!(\w+)\][+-]?\s*(.*)$/);
    if (!m) { out.push(line); continue; }
    const body = [];
    while (i + 1 < lines.length && /^>/.test(lines[i + 1])) body.push(lines[++i].replace(/^>\s?/, ""));
    const type = m[1].toLowerCase();
    const title = m[2] || type[0].toUpperCase() + type.slice(1);
    // The title sits inside raw HTML, so render its math here rather than via placeholders.
    const titleHtml = title.split(/(\$[^$\n]+\$)/g).map((s, k) =>
      k % 2 ? katex.renderToString(s.slice(1, -1), { throwOnError: false, strict: "ignore" }) : md.renderInline(s)).join("");
    out.push("", `<div class="callout callout-${type}"><div class="callout-title">${titleHtml}</div>`, "",
      convertCallouts(body.join("\n")), "", "</div>", "");
  }
  return out.join("\n");
}

function chapterHtml(rel, blocks, mathErrors) {
  const abs = path.join(ROOT, rel);
  const dir = path.dirname(abs);
  const id = docId(rel);
  let text = fs.readFileSync(abs, "utf8").replace(/\r\n/g, "\n");

  // Desmos blocks -> cached screenshots.
  text = text.replace(/```desmos-graph\n([\s\S]*?)```/g, (_, body) => {
    const h = desmosHash(body);
    blocks.set(h, body);
    if (!blockFiles.has(h)) blockFiles.set(h, new Set());
    blockFiles.get(h).add(rel);
    return `<img class="desmos" src="${pathToFileURL(path.join(CACHE, h + ".png")).href}" alt="Desmos graph">`;
  });
  // Drop the in-file navigation lines (prev/next); the PDF has its own outline.
  text = text.replace(/^\[← [^\n]*\]\([^\n]*\)[^\n]*$/gm, (l) => (/(Previous|Table of Contents|README|Next)/.test(l) ? "" : l));
  text = convertCallouts(text);
  const store = [];
  text = extractMath(text, store);

  const env = {};
  const tokens = md.parse(text, env);
  // Heading ids: "<docid>--<slug>"; the first H1 is the document anchor itself.
  let first = true;
  for (let i = 0; i < tokens.length; i++) {
    const t = tokens[i];
    if (t.type === "heading_open") {
      const inline = tokens[i + 1];
      const plain = inline.children.filter((c) => c.type === "text" || c.type === "code_inline").map((c) => c.content).join("");
      t.attrSet("id", first && t.tag === "h1" ? id : `${id}--${slugify(plain.replace(/KTXPH\d+ENDPH/g, ""))}`);
      if (first && t.tag === "h1") first = false;
    }
    // Rewrite links and image sources.
    const fix = (c) => {
      if (c.type === "link_open") {
        const href = c.attrGet("href") || "";
        if (/^(https?:|mailto:)/.test(href)) return;
        const [p, frag] = href.split("#");
        if (!p) { c.attrSet("href", `#${id}--${frag}`); return; }
        if (!p.endsWith(".md")) { c.attrSet("href", pathToFileURL(path.resolve(dir, p)).href); return; }
        const target = path.relative(ROOT, path.resolve(dir, decodeURIComponent(p))).split(path.sep).join("/");
        c.attrSet("href", frag ? `#${docId(target)}--${frag}` : `#${docId(target)}`);
      } else if (c.type === "image") {
        const src = c.attrGet("src");
        if (src && !/^(https?:|file:|data:)/.test(src)) c.attrSet("src", pathToFileURL(path.resolve(dir, decodeURIComponent(src))).href);
      }
    };
    if (t.children) t.children.forEach(fix);
  }
  let html = md.renderer.render(tokens, md.options, env);
  html = renderMath(html, store, mathErrors, rel);
  const kind = rel.startsWith("tracks/") && rel.endsWith("README.md") ? "track-index" : "chapter";
  return `<section class="${kind}" data-src="${rel}">${html}</section>`;
}

// ---------------------------------------------------------------- assemble

function titleOf(rel) {
  const first = fs.readFileSync(path.join(ROOT, rel), "utf8").split("\n")[0];
  return first.replace(/^#\s*/, "");
}

function buildHtml(items, blocks, mathErrors) {
  const body = [];
  const toc = [];
  for (const it of items) {
    if (it.divider) {
      const id = "div-" + slugify(it.divider);
      body.push(`<section class="divider${it.minor ? " minor" : ""}"><h1 id="${id}">${escapeHtml(it.divider)}</h1>${it.blurb ? `<p>${escapeHtml(it.blurb)}</p>` : ""}</section>`);
      toc.push(`<li class="${it.minor ? "toc-part" : "toc-level"}"><a href="#${id}">${escapeHtml(it.divider)}</a></li>`);
    } else {
      body.push(chapterHtml(it.file, blocks, mathErrors));
      toc.push(`<li class="toc-ch"><a href="#${docId(it.file)}">${md.renderInline(titleOf(it.file))}</a></li>`);
    }
  }
  const katexCss = pathToFileURL(path.join(HERE, "node_modules", "katex", "dist", "katex.min.css")).href;
  const style = fs.readFileSync(path.join(HERE, "book.css"), "utf8");
  const cover = pathToFileURL(path.join(ROOT, "images", "universe-map.svg")).href;
  return `<!doctype html><html lang="en"><head><meta charset="utf-8"><title>${BOOK.title}</title>
<link rel="stylesheet" href="${katexCss}"><style>${style}</style></head><body>
<section class="cover"><h1 class="cover-title">${BOOK.title}</h1><p class="cover-sub">${BOOK.subtitle}</p>
${fs.existsSync(path.join(ROOT, "images", "universe-map.svg")) ? `<img class="cover-map" src="${cover}" alt="">` : ""}
<p class="cover-note">Levels 1–3 · Core, specialization tracks, research frontier · Appendix of proofs</p></section>
<section class="toc"><h1 id="contents">Contents</h1><ul>${toc.join("\n")}</ul></section>
${body.join("\n")}
</body></html>`;
}

/** Chrome writes only the title; fill in the rest of the document properties from book.json. */
async function setMetadata(file) {
  const pdf = await PDFDocument.load(fs.readFileSync(file), { updateMetadata: false });
  pdf.setTitle(ONLY ? `${BOOK.title} (${ONLY})` : BOOK.title, { showInWindowTitleBar: true });
  pdf.setAuthor(BOOK.author);
  pdf.setSubject(`${BOOK.subtitle}. ${BOOK.description}`);
  pdf.setKeywords([BOOK.keywords.join(", ")]);
  pdf.setCreator("pdf/build.mjs (Markdown, KaTeX, Desmos, headless Chrome)");
  pdf.setProducer("Skia/PDF via Chrome; metadata by pdf-lib");
  pdf.setLanguage("en");
  const now = new Date();
  pdf.setCreationDate(now);
  pdf.setModificationDate(now);
  fs.writeFileSync(file, await pdf.save({ useObjectStreams: true }));
}

async function launch() {
  for (const channel of ["chrome", "msedge"]) {
    try { return await chromium.launch({ channel }); } catch { /* try next */ }
  }
  throw new Error("Need Google Chrome or Microsoft Edge installed (playwright-core uses the system browser).");
}

async function main() {
  const items = readingOrder();
  const blocks = new Map();
  const mathErrors = [];
  let html = buildHtml(items, blocks, mathErrors);
  const browser = await launch();
  try {
    await renderDesmos(blocks, browser);
    reportDesmosErrors(blocks);
    if (args.includes("--desmos-only")) return;
    fs.mkdirSync(path.dirname(OUT_HTML), { recursive: true });
    fs.writeFileSync(OUT_HTML, html);
    console.log(`html: ${items.filter((i) => i.file).length} documents, ${(html.length / 1e6).toFixed(1)} MB; ${mathErrors.length} KaTeX errors`);
    fs.writeFileSync(path.join(HERE, "cache", "katex-errors.txt"), mathErrors.join("\n") + (mathErrors.length ? "\n" : ""));
    html = null;
    const page = await browser.newPage();
    page.setDefaultTimeout(0);
    await page.goto(pathToFileURL(OUT_HTML).href, { waitUntil: "load", timeout: 0 });
    await page.evaluate(() => document.fonts.ready);
    console.log("printing PDF ...");
    await page.pdf({
      path: OUT_PDF, format: "A4", printBackground: true, outline: true, tagged: true,
      margin: { top: "18mm", bottom: "18mm", left: "17mm", right: "17mm" },
      displayHeaderFooter: true,
      headerTemplate: `<div style="font-size:7.5px;color:#9ca3af;width:100%;text-align:right;padding-right:17mm">${BOOK.title}</div>`,
      footerTemplate: `<div style="font-size:8px;color:#6b7280;width:100%;text-align:center"><span class="pageNumber"></span></div>`,
    });
    await setMetadata(OUT_PDF);
    console.log(`pdf: ${OUT_PDF} (${(fs.statSync(OUT_PDF).size / 1e6).toFixed(1)} MB)`);
  } finally {
    await browser.close();
  }
}

main().catch((e) => { console.error(e); process.exit(1); });
