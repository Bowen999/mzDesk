import { Link, useInRouterContext } from "react-router-dom";
import { BASE } from "../api";
import { Icon } from "./Icon";
import { RichText } from "./Help";

/** The one-line description under a title: "what I have → what I get". The arrow splits it; the result side is the stronger one. */
export function Flow({ text }: { text: string }) {
  const at = text.indexOf("→");
  if (at < 0) return <p className="one-liner"><RichText text={text} /></p>;
  return (
    <p className="one-liner flow">
      <span className="flow-in"><RichText text={text.slice(0, at).trim()} /></span>
      <span className="flow-to">
        <span className="flow-arrow" aria-hidden="true">→</span>
        <span className="sr-only"> gives </span>
        <span className="flow-out"><RichText text={text.slice(at + 1).trim()} /></span>
      </span>
    </p>
  );
}

/** Title block of every page (module or workflow): icon, title, the "have → get" line. No kicker, no subtitle. */
export function PageHead({ icon, title, flow, preview = false }: { icon: string; title: string; flow: string; preview?: boolean }) {
  return (
    <div className="module-heading">
      <span className="module-emblem"><Icon name={icon} /></span>
      <div className="module-heading-copy">
        <h1>{title}{preview && <span className="badge tone-neutral head-badge">Preview</span>}</h1>
        <Flow text={flow} />
      </div>
    </div>
  );
}

/** Foot of every page: what it is built on, and the link to its long notes. `link` is off in the Session report (a file has no app to open). */
export function PageFoot({ id, basedOn, link = true }: { id: string; basedOn?: string[]; link?: boolean }) {
  const inRouter = useInRouterContext();
  if (!basedOn?.length && !link) return null;
  return (
    <footer className="page-foot">
      {!!basedOn?.length && <span>Based on {basedOn.join(", ")}</span>}
      {link && (inRouter ? <Link to={`/docs/${id}`}>Docs</Link> : <a href={`${BASE}docs/${id}`}>Docs</a>)}
    </footer>
  );
}
