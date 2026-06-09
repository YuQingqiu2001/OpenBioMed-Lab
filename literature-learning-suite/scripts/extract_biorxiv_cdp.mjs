#!/usr/bin/env node
/**
 * bioRxiv/medRxiv browser full-text fetcher via Chrome DevTools Protocol.
 *
 * 设计原则：
 * - 使用用户可见的真实 Chrome + 用户手动通过 Cloudflare/Turnstile；不做验证码破解或指纹伪装。
 * - 自动完成“打开页面 -> 等待渲染 -> 提取正文 -> 写入 fulltext_cache”。
 * - 适合少量精读/半自动下载；批量监控仍建议使用 download_biorxiv_api.py。
 *
 * Usage:
 *   node scripts/extract_biorxiv_cdp.mjs --doi 10.64898/2026.05.31.727600 --port 9223
 *   node scripts/extract_biorxiv_cdp.mjs --server medrxiv --doi 10.1101/2026.01.01.123456 --port 9223
 *   node scripts/extract_biorxiv_cdp.mjs --url https://www.biorxiv.org/content/10.64898/2026.05.31.727600v1.full --port 9223
 *   node scripts/extract_biorxiv_cdp.mjs --batch ./doi_list.txt --port 9223
 */

import { mkdirSync, readFileSync, writeFileSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

const DEFAULT_CDP = 'http://127.0.0.1:9223';
const SCRIPT_DIR = dirname(fileURLToPath(import.meta.url));
const DEFAULT_ROOT = process.env.LITERATURE_KG_ROOT || join(SCRIPT_DIR, '..', 'literature-workspace');
const DEFAULT_OUT_DIR = join(DEFAULT_ROOT, 'fulltext_cache');

function usage(exitCode = 0) {
  console.log(`bioRxiv/medRxiv browser full-text fetcher

Options:
  --doi DOI              Single DOI, e.g. 10.64898/2026.05.31.727600
  --url URL              Direct article URL
  --batch FILE           Text file: one DOI or URL per line (# comments allowed)
  --server NAME          biorxiv or medrxiv (default: biorxiv)
  --version N            Version number when DOI API cannot resolve it (default: 1)
  --port N               Chrome CDP port (default: 9223)
  --cdp URL              Full CDP URL (default: http://127.0.0.1:9223)
  --out-dir DIR          Output directory (default: LITERATURE_KG_ROOT/fulltext_cache)
  --wait-human SEC       Seconds to wait for manual Cloudflare verification (default: 180)
  --delay SEC            Delay between batch items (default: 8)
  --keep-open            Do not close created tab after extraction
  --debug-html           Save full rendered HTML beside the text file
  --help                 Show help

Important:
  This script does not bypass Cloudflare or solve CAPTCHAs. If a security check appears,
  complete it manually in the visible Chrome window, then the script continues.
`);
  process.exit(exitCode);
}

function parseArgs(argv) {
  const opts = {
    server: 'biorxiv',
    version: null,
    cdp: null,
    port: '9223',
    outDir: DEFAULT_OUT_DIR,
    waitHuman: 180,
    delay: 8,
    keepOpen: false,
    debugHtml: false,
  };
  const positional = [];
  for (let i = 2; i < argv.length; i++) {
    const a = argv[i];
    if (a === '--help' || a === '-h') usage(0);
    if (a === '--keep-open') { opts.keepOpen = true; continue; }
    if (a === '--debug-html') { opts.debugHtml = true; continue; }
    if (a.startsWith('--')) {
      const key = a.slice(2).replace(/-([a-z])/g, (_, c) => c.toUpperCase());
      const val = argv[++i];
      if (val == null) throw new Error(`Missing value for ${a}`);
      opts[key] = val;
    } else {
      positional.push(a);
    }
  }
  if (!opts.doi && !opts.url && !opts.batch && positional.length) {
    const p = positional[0];
    if (/^https?:\/\//i.test(p)) opts.url = p;
    else opts.doi = p;
  }
  opts.server = String(opts.server || 'biorxiv').toLowerCase();
  if (!['biorxiv', 'medrxiv'].includes(opts.server)) {
    throw new Error('--server must be biorxiv or medrxiv');
  }
  opts.waitHuman = Number(opts.waitHuman ?? 180);
  opts.delay = Number(opts.delay ?? 8);
  opts.cdp = opts.cdp || `http://127.0.0.1:${opts.port || '9223'}`;
  return opts;
}

function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }
function safeName(s) { return String(s).replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, ''); }
function stripKnownPrefix(s) { return String(s || '').replace(/^(BIORXIV|MEDRXIV):/i, '').trim(); }

function parseDoiAndVersion(input) {
  let doi = stripKnownPrefix(input).trim();
  doi = doi.replace(/^doi:\s*/i, '').trim();
  const mUrl = doi.match(/\/content\/(10\.[^\s?#]+?)(?:\.full|\.full\.pdf|\.source\.xml)?(?:[?#].*)?$/i);
  if (mUrl) doi = decodeURIComponent(mUrl[1]);
  const m = doi.match(/^(10\.\d+\/.+?)(?:v(\d+))?$/i);
  if (!m) return { doi, version: null };
  return { doi: m[1], version: m[2] || null };
}

async function cdpJson(url, options = {}) {
  const r = await fetch(url, options);
  if (!r.ok) throw new Error(`${url} -> HTTP ${r.status}`);
  return await r.json();
}

async function getBrowserVersion(cdp) {
  return await cdpJson(`${cdp}/json/version`);
}

async function newTab(cdp) {
  const r = await fetch(`${cdp}/json/new`, { method: 'PUT' });
  if (!r.ok) throw new Error(`${cdp}/json/new -> HTTP ${r.status}`);
  return await r.json();
}

async function closeTab(cdp, id) {
  try { await fetch(`${cdp}/json/close/${id}`); } catch (_) {}
}

function connectWs(url) {
  return new Promise((resolve, reject) => {
    const ws = new WebSocket(url);
    ws.onopen = () => resolve(ws);
    ws.onerror = () => reject(new Error('WebSocket connection failed'));
  });
}

function makeSender(ws) {
  let id = 0;
  return function send(method, params = {}, timeoutMs = 30000) {
    const msgId = ++id;
    return new Promise((resolve, reject) => {
      const t = setTimeout(() => {
        ws.removeEventListener('message', onMessage);
        reject(new Error(`${method} timeout`));
      }, timeoutMs);
      function onMessage(event) {
        const msg = JSON.parse(event.data);
        if (msg.id !== msgId) return;
        clearTimeout(t);
        ws.removeEventListener('message', onMessage);
        if (msg.error) reject(new Error(`${method}: ${msg.error.message}`));
        else resolve(msg.result ?? msg);
      }
      ws.addEventListener('message', onMessage);
      ws.send(JSON.stringify({ id: msgId, method, params }));
    });
  };
}

async function waitForLoad(ws, timeoutMs = 45000) {
  return new Promise(resolve => {
    const t = setTimeout(() => {
      ws.removeEventListener('message', onMessage);
      resolve(false);
    }, timeoutMs);
    function onMessage(event) {
      const msg = JSON.parse(event.data);
      if (msg.method === 'Page.loadEventFired') {
        clearTimeout(t);
        ws.removeEventListener('message', onMessage);
        resolve(true);
      }
    }
    ws.addEventListener('message', onMessage);
  });
}

async function evaluateJson(send, expression, timeoutMs = 30000) {
  const res = await send('Runtime.evaluate', {
    expression,
    awaitPromise: true,
    returnByValue: true,
  }, timeoutMs);
  if (res?.exceptionDetails) {
    const detail = JSON.stringify(res.exceptionDetails, null, 2);
    const msg = res.exceptionDetails.exception?.description || res.exceptionDetails.text || detail;
    throw new Error(`Runtime.evaluate failed: ${msg}\n${detail}`);
  }
  const remote = res?.result;
  if (!remote) throw new Error(`Runtime.evaluate returned no result: ${JSON.stringify(res).slice(0, 500)}`);
  const value = remote.value;
  if (typeof value === 'string') return JSON.parse(value);
  if (value === undefined && remote.type === 'undefined') {
    throw new Error(`Runtime.evaluate returned undefined for expression: ${expression.slice(0, 160)}`);
  }
  return value;
}

const pageStatusExpr = String.raw`JSON.stringify((function(){
  const body = document.body ? document.body.innerText : '';
  const head = body.slice(0, 1200);
  const title = document.title || '';
  const hay = (title + '\n' + head).toLowerCase();
  const isCloudflare = /just a moment|请稍候|安全验证|verify you are human|checking if the site connection is secure|turnstile|cloudflare/.test(hay);
  const articleSelectors = ['.article.fulltext-view', '#fulltext-view', '.article-body', '.section.abstract', 'article', '#content-block', 'main'];
  const selector = articleSelectors.find(s => document.querySelector(s));
  const hasArticleWords = /\b(Abstract|Introduction|Methods|Results|Discussion|References)\b/i.test(body);
  return {
    url: location.href,
    title,
    isCloudflare,
    selector: selector || '',
    textLength: body.length,
    hasArticle: Boolean(selector && body.length > 1500 && hasArticleWords && !isCloudflare),
    preview: head
  };
})())`;

const extractExpr = String.raw`JSON.stringify((function(){
  function text(sel) {
    const el = document.querySelector(sel);
    return el ? (el.innerText || el.textContent || '').trim() : '';
  }
  function firstText(selectors) {
    for (const s of selectors) {
      const t = text(s);
      if (t && t.length > 20) return { selector: s, text: t };
    }
    return { selector: 'body', text: document.body ? document.body.innerText : '' };
  }
  const title = firstText(['h1.highwire-cite-title', 'h1.citation__title', 'h1', '.highwire-cite-title']).text;
  const authors = firstText(['.highwire-citation-authors', '.authors', '.contrib-group', '.article-authors']).text;
  const doi = (document.body ? (document.body.innerText.match(/10\.\d+\/[-._;()/:A-Z0-9]+/i) || [''])[0] : '');
  const candidates = [
    '.article.fulltext-view',
    '#fulltext-view',
    '.article-body',
    '#content-block .article',
    '#content-block',
    'article',
    'main',
    'body'
  ];
  const picked = firstText(candidates);
  let raw = picked.text || '';
  const lines = raw.split(/\n+/).map(x => x.trim()).filter(Boolean);
  const drop = [
    /^Skip to main content$/i, /^Log in$/i, /^Sign in$/i, /^Search$/i, /^Menu$/i,
    /^Download PDF$/i, /^PDF$/i, /^Email alerts$/i, /^Share$/i, /^Article usage$/i,
    /^Copyright/i, /^Cold Spring Harbor Laboratory/i, /^bioRxiv/i, /^medRxiv/i
  ];
  const filtered = [];
  for (const line of lines) {
    if (drop.some(r => r.test(line))) continue;
    if (filtered.length && filtered[filtered.length - 1] === line) continue;
    filtered.push(line);
  }
  const links = Array.from(document.querySelectorAll('a[href]')).map(a => ({
    text: (a.innerText || '').trim(),
    href: new URL(a.getAttribute('href'), location.href).href
  })).filter(x => /\.full\.pdf|\.source\.xml|download|pdf|xml/i.test(x.href + ' ' + x.text)).slice(0, 30);
  return {
    url: location.href,
    title,
    authors,
    doi,
    selector: picked.selector,
    text: filtered.join('\n'),
    rawLength: raw.length,
    cleanedLength: filtered.join('\n').length,
    links,
    html: document.documentElement.outerHTML
  };
})())`;

async function fetchApiRecord(server, doi) {
  try {
    const url = `https://api.biorxiv.org/details/${server}/${doi}/na/json`;
    const payload = await cdpJson(url);
    const rec = payload?.collection?.[0] || null;
    if (rec) console.log(`[api] ${rec.title || ''} | version=${rec.version || '?'} | date=${rec.date || '?'}`);
    return rec;
  } catch (e) {
    console.log(`[api] DOI lookup failed: ${e.message}`);
    return null;
  }
}

function buildUrl(server, doi, version, directUrl) {
  if (directUrl) return directUrl;
  return `https://www.${server}.org/content/${doi}v${version}.full`;
}

async function processOne(input, opts, browserInfo) {
  let server = opts.server;
  let url = opts.url && !opts.batch ? opts.url : null;
  let doi = opts.doi || input || '';
  let version = opts.version || null;

  if (url || /^https?:\/\//i.test(doi)) {
    url = url || doi;
    server = /medrxiv\.org/i.test(url) ? 'medrxiv' : 'biorxiv';
    const parsed = parseDoiAndVersion(url);
    doi = parsed.doi;
    version = version || parsed.version;
  } else {
    const parsed = parseDoiAndVersion(doi);
    doi = parsed.doi;
    version = version || parsed.version;
  }
  if (!doi && !url) throw new Error('No DOI or URL provided');

  const rec = doi ? await fetchApiRecord(server, doi) : null;
  version = version || rec?.version || '1';
  url = buildUrl(server, doi, version, url);

  const prefix = server.toUpperCase();
  // Cache convention: keep DOI-based fulltext cache stable across retries.
  // Version is recorded in metadata/header; do not suffix _v1 by default because
  // gen_edges.py and other readers look for normalized paper IDs like
  // BIORXIV_10.64898_2026.05.31.727600.txt.
  const baseName = doi ? safeName(`${prefix}_${doi}`) : safeName(`${prefix}_${url}_v${version}`);
  const outPath = join(opts.outDir, `${baseName}.txt`);
  const jsonPath = join(opts.outDir, `${baseName}.metadata.json`);
  const htmlPath = join(opts.outDir, `${baseName}.html`);
  mkdirSync(dirname(outPath), { recursive: true });

  console.log(`\n=== Fetch ${prefix}:${doi || url} ===`);
  console.log(`[target] ${url}`);

  const tab = await newTab(opts.cdp);
  const ws = await connectWs(tab.webSocketDebuggerUrl);
  const send = makeSender(ws);
  try {
    await send('Page.enable');
    await send('Runtime.enable');
    await send('Network.enable');

    await send('Page.navigate', { url }, 10000);
    await waitForLoad(ws, 45000);
    await sleep(3000);

    const deadline = Date.now() + opts.waitHuman * 1000;
    let warned = false;
    let status = null;
    while (Date.now() < deadline) {
      status = await evaluateJson(send, pageStatusExpr, 10000);
      console.log(`[status] title=${JSON.stringify(status.title).slice(0, 90)} len=${status.textLength} cf=${status.isCloudflare} selector=${status.selector || '-'} url=${status.url}`);
      if (status.hasArticle) break;
      if (status.isCloudflare && !warned) {
        console.log('\n[manual action needed] Chrome 中出现 Cloudflare/安全验证。请在打开的可见 Chrome 窗口里手动完成验证；脚本会继续等待。\n');
        warned = true;
      }
      if (!status.isCloudflare && status.textLength > 3000 && /Abstract|Introduction|References/i.test(status.preview)) break;
      await sleep(3000);
    }

    const extracted = await evaluateJson(send, extractExpr, 30000);
    const text = extracted.text || '';
    const ok = text.length > 2000 && !/Just a moment|正在进行安全验证|Verify you are human/i.test(text.slice(0, 1000));

    const header = [
      `# Downloaded: ${new Date().toISOString()}`,
      `# Source: browser-cdp`,
      `# URL: ${extracted.url || url}`,
      `# Server: ${server}`,
      `# DOI: ${doi || extracted.doi || ''}`,
      `# Version: ${version}`,
      `# Title: ${(extracted.title || rec?.title || '').replace(/\s+/g, ' ').trim()}`,
      `# Browser: ${browserInfo.Browser || ''}`,
      `# Note: Extracted from a user-visible Chrome session; no CAPTCHA bypass or fingerprint spoofing.`,
      '',
    ].join('\n');

    const textPath = ok ? outPath : outPath.replace(/\.txt$/i, '.failed.txt');
    writeFileSync(textPath, header + text + '\n', 'utf-8');
    writeFileSync(jsonPath, JSON.stringify({
      ok,
      server,
      doi,
      version,
      requestedUrl: url,
      finalUrl: extracted.url,
      title: extracted.title || rec?.title || '',
      authors: extracted.authors || rec?.authors || '',
      apiRecord: rec,
      selector: extracted.selector,
      rawLength: extracted.rawLength,
      cleanedLength: extracted.cleanedLength,
      links: extracted.links,
      outPath: textPath,
      intendedOutPath: outPath,
      extractedAt: new Date().toISOString(),
    }, null, 2), 'utf-8');
    if (opts.debugHtml) writeFileSync(htmlPath, extracted.html || '', 'utf-8');

    console.log(`[saved] ${textPath}`);
    console.log(`[meta]  ${jsonPath}`);
    if (opts.debugHtml) console.log(`[html]  ${htmlPath}`);
    console.log(`[result] ok=${ok} cleanedLength=${text.length} selector=${extracted.selector}`);
    console.log('[preview]');
    console.log(text.slice(0, 1200));

    return { ok, outPath: textPath, jsonPath, length: text.length, title: extracted.title || rec?.title || '' };
  } finally {
    try { ws.close(); } catch (_) {}
    if (!opts.keepOpen) await closeTab(opts.cdp, tab.id);
  }
}

async function main() {
  const opts = parseArgs(process.argv);
  if (!opts.doi && !opts.url && !opts.batch) usage(1);

  const browserInfo = await getBrowserVersion(opts.cdp);
  console.log(`[cdp] ${opts.cdp}`);
  console.log(`[browser] ${browserInfo.Browser || '?'} | UA=${browserInfo['User-Agent'] || '?'}`);
  if (/HeadlessChrome/i.test(browserInfo['User-Agent'] || '')) {
    console.log('\n[warning] 当前 CDP 是 HeadlessChrome。bioRxiv 全文/Cloudflare 对 headless 很不友好；建议运行 biorxiv_chrome_cdp_launcher.bat 打开可见 Chrome，再使用 --port 9223。\n');
  }

  if (opts.batch) {
    const lines = readFileSync(opts.batch, 'utf-8')
      .split(/\r?\n/)
      .map(x => x.trim())
      .filter(x => x && !x.startsWith('#'));
    console.log(`[batch] ${lines.length} items from ${opts.batch}`);
    const results = [];
    for (let i = 0; i < lines.length; i++) {
      console.log(`\n[batch] ${i + 1}/${lines.length}`);
      try {
        results.push(await processOne(lines[i], { ...opts, doi: null, url: null }, browserInfo));
      } catch (e) {
        console.error(`[error] ${lines[i]} -> ${e.stack || e.message}`);
        results.push({ ok: false, input: lines[i], error: e.message });
      }
      if (i + 1 < lines.length) await sleep(opts.delay * 1000);
    }
    const summaryPath = join(opts.outDir, `BIORXIV_BATCH_${new Date().toISOString().replace(/[:.]/g, '-')}.json`);
    writeFileSync(summaryPath, JSON.stringify(results, null, 2), 'utf-8');
    console.log(`\n[batch summary] ${summaryPath}`);
    const okCount = results.filter(x => x.ok).length;
    console.log(`[batch result] ok=${okCount}/${results.length}`);
    return okCount === results.length ? 0 : 2;
  }

  const result = await processOne(opts.doi || opts.url, opts, browserInfo);
  return result.ok ? 0 : 2;
}

main().then(code => process.exit(code)).catch(e => {
  console.error(`[fatal] ${e.stack || e.message}`);
  process.exit(1);
});
