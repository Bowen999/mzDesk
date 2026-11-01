import { useEffect, useState } from "react";
import { Link, useParams } from "react-router-dom";
import { getDocs } from "../api";
import { Markdown } from "../components/Markdown";
import type { Spec } from "../types";

/** /docs/<page id>: the long notes of a module or workflow (mzdesk/docs/<id>.md). The pages keep one short line; the rest lives here. */
export function DocsPage({ spec }: { spec: Spec }) {
  const { id = "" } = useParams();
  const module = spec.modules.find(m => m.id === id), pipeline = spec.pipelines.find(p => p.id === id);
  const page = module ?? pipeline;
  const [text, setText] = useState<string | null>(null);
  const [failed, setFailed] = useState(false);
  useEffect(() => {
    let live = true;
    setText(null); setFailed(false);
    getDocs(id).then(d => { if (live) setText(d.markdown); }).catch(() => { if (live) setFailed(true); });
    return () => { live = false; };
  }, [id]);
  if (!page) return <p className="page">Unknown page.</p>;
  const back = module ? `/modules/${id}` : `/pipelines/${id}`;
  return (
    <div className="page doc-page">
      <Link className="doc-back" to={back}>← {page.title}</Link>
      {text === null && !failed && <div className="doc-loading" aria-hidden="true"><i /><i /><i /></div>}
      {failed && <p className="msg">No notes for this page yet.</p>}
      {text !== null && <Markdown text={text} />}
    </div>
  );
}
