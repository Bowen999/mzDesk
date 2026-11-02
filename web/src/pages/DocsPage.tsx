import { useEffect, useMemo, useState } from "react";
import { Link, Navigate, NavLink, useLocation } from "react-router-dom";
import { loadPage } from "../docs/load";
import { docsOutline, flatPages, hasTitle, headingsOf, isFolder, type DocItem, type DocSection } from "../docs/outline";
import { Logo } from "../components/Logo";
import { Markdown } from "../components/Markdown";
import { ThemeSwitch } from "../components/Sidebar";
import type { Spec } from "../types";
import "../docs/docs.css";

/** /docs and /docs/<page id>: the Docs, a page of its own (top bar, contents on the left, the page, "On this page" on the right).
 *  English only. Every page shows its text from /api/docs/<id> (mzdesk/docs/<id>.md).
 *  `backTo` is the app page the reader came from. */
export function DocsLayout({ spec, backTo }: { spec: Spec; backTo: string }) {
  const { pathname } = useLocation();
  const id = pathname.split("/")[2] ?? "";
  const outline = useMemo(() => docsOutline(spec), [spec]);
  const pages = useMemo(() => flatPages(outline), [outline]);
  const [navOpen, setNavOpen] = useState(false);
  useEffect(() => setNavOpen(false), [pathname]);
  useEffect(() => {
    if (!navOpen) return;
    const key = (e: KeyboardEvent) => { if (e.key === "Escape") setNavOpen(false); };
    document.addEventListener("keydown", key);
    return () => document.removeEventListener("keydown", key);
  }, [navOpen]);

  if (!id) return <Navigate to={`/docs/${pages[0].id}`} replace />;
  const at = pages.findIndex(p => p.id === id);
  const item = pages[at];
  const section = outline.find(s => s.entries.some(e => (isFolder(e) ? e.items.some(i => i.id === id) : e.id === id)));

  return (
    <div className={"docs" + (navOpen ? " nav-open" : "")}>
      <header className="docs-top">
        <button type="button" className="docs-menu" aria-label="Open contents" aria-expanded={navOpen} aria-controls="docs-nav" onClick={() => setNavOpen(true)}>
          <span aria-hidden="true" /><span aria-hidden="true" /><span aria-hidden="true" />
        </button>
        <Link to={`/docs/${pages[0].id}`} className="docs-brand"><Logo /><span className="docs-brand-name">mzDesk</span><span className="docs-brand-tag">Docs</span></Link>
        <div className="docs-top-actions">
          <Link to={backTo} className="docs-back">← Back to mzDesk</Link>
          <ThemeSwitch compact />
        </div>
      </header>
      <div className="docs-backdrop" aria-hidden="true" onClick={() => setNavOpen(false)} />
      <DocsNav outline={outline} activeId={id} />
      <main className="docs-main">
        {item ? <DocsArticle key={item.id} item={item} section={section?.title} prev={pages[at - 1]} next={pages[at + 1]} />
          : <article className="docs-article"><h1>Page not found</h1><p className="docs-empty">There is no Docs page called “{id}”.</p></article>}
      </main>
    </div>
  );
}

/** The contents on the left: one block per section; module subgroups fold (the one holding the open page starts open). */
function DocsNav({ outline, activeId }: { outline: DocSection[]; activeId: string }) {
  const [open, setOpen] = useState<Record<string, boolean>>({});
  const link = (i: DocItem, sub = false) => (
    <NavLink key={i.id} to={`/docs/${i.id}`} className={({ isActive }) => "docs-link" + (sub ? " sub" : "") + (isActive ? " active" : "")}>{i.title}</NavLink>
  );
  return (
    <nav className="docs-nav" id="docs-nav" aria-label="Docs contents">
      {outline.map(s => (
        <div key={s.title} className="docs-nav-section">
          <div className="docs-nav-title">{s.title}</div>
          {s.entries.map(e => {
            if (!isFolder(e)) return link(e);
            const hasActive = e.items.some(i => i.id === activeId);
            const isOpen = open[e.id] ?? hasActive;
            return (
              <div key={e.id} className="docs-folder">
                <button type="button" className={"docs-link docs-folder-title" + (hasActive ? " has-active" : "")} aria-expanded={isOpen}
                  onClick={() => setOpen({ ...open, [e.id]: !isOpen })}>
                  <span>{e.title}</span><span className={"docs-caret" + (isOpen ? " open" : "")} aria-hidden="true">›</span>
                </button>
                {isOpen && <div className="docs-folder-items">{e.items.map(i => link(i, true))}</div>}
              </div>
            );
          })}
        </div>
      ))}
    </nav>
  );
}

/** One page: title, the text (/api/docs/<id>), a link into the app, Previous / Next, and "On this page". */
function DocsArticle({ item, section, prev, next }: { item: DocItem; section?: string; prev?: DocItem; next?: DocItem }) {
  const { hash, key } = useLocation();
  const [text, setText] = useState<string | null>(null);
  const [problem, setProblem] = useState<"missing" | "failed" | null>(null);
  useEffect(() => {
    let live = true;
    loadPage(item.id).then(t => { if (!live) return; if (t === null) setProblem("missing"); else setText(t); }).catch(() => { if (live) setProblem("failed"); });
    return () => { live = false; };
  }, [item.id]);
  // a new page starts at the top; a link with #anchor scrolls to its heading once the text is there
  useEffect(() => {
    if (!hash) { if (text === null && !problem) window.scrollTo(0, 0); return; }
    if (text === null) return;
    let id = hash.slice(1);
    try { id = decodeURIComponent(id); } catch { /* keep the raw id */ }
    document.getElementById(id)?.scrollIntoView({ block: "start" });
  }, [hash, key, text, problem]);
  const headings = useMemo(() => (text ? headingsOf(text) : []), [text]);
  const appPath = item.kind === "module" ? `/modules/${item.id}` : item.kind === "workflow" ? `/pipelines/${item.id}` : null;

  return (
    <div className="docs-page">
      <article className="docs-article">
        {section && <div className="docs-eyebrow">{section}</div>}
        {(text === null || !hasTitle(text)) && <h1>{item.title}</h1>}
        {text === null && !problem && <div className="doc-loading" aria-hidden="true"><i /><i /><i /></div>}
        {text !== null && <Markdown text={text} />}
        {problem === "missing" && <p className="docs-empty">This page is not written yet.</p>}
        {problem === "failed" && <p className="docs-empty">This page could not be loaded. Reload to try again.</p>}
        {appPath && <p className="docs-open"><Link to={appPath}>Open {item.title} in mzDesk →</Link></p>}
        <nav className="docs-pager" aria-label="Previous and next page">
          {prev ? <Link to={`/docs/${prev.id}`} className="docs-pager-link"><span>Previous</span><b>{prev.title}</b></Link> : <span />}
          {next ? <Link to={`/docs/${next.id}`} className="docs-pager-link next"><span>Next</span><b>{next.title}</b></Link> : <span />}
        </nav>
      </article>
      <OnThisPage headings={headings} />
    </div>
  );
}

/** "On this page": the ## / ### headings; the one last scrolled past is marked. */
function OnThisPage({ headings }: { headings: { id: string; text: string; level: 2 | 3 }[] }) {
  const [active, setActive] = useState("");
  useEffect(() => {
    if (!headings.length) return;
    const spy = () => {
      let cur = headings[0].id;
      for (const h of headings) {
        const el = document.getElementById(h.id);
        if (el && el.getBoundingClientRect().top < 110) cur = h.id;
      }
      setActive(cur);
    };
    spy();
    window.addEventListener("scroll", spy, { passive: true });
    return () => window.removeEventListener("scroll", spy);
  }, [headings]);
  return (
    <aside className="docs-toc" aria-label="On this page">
      {!!headings.length && <>
        <div className="docs-toc-title">On this page</div>
        {headings.map(h => (
          <a key={h.id} href={`#${h.id}`} className={"docs-toc-link" + (h.level === 3 ? " sub" : "") + (active === h.id ? " active" : "")}
            onClick={e => { e.preventDefault(); document.getElementById(h.id)?.scrollIntoView({ behavior: "smooth", block: "start" }); }}>{h.text}</a>
        ))}
      </>}
    </aside>
  );
}
