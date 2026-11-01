import type { ReactNode } from "react";
import { RichText } from "./Help";

/** The small markdown subset of the Docs files (mzdesk/docs/*.md): #, ##, ### headings, "- " bullets, paragraphs, **bold**, `code`
 *  and the <sub> <sup> <i> <b> tags that RichText draws. Nothing else is interpreted. */
function inline(text: string): ReactNode {
  return text.split(/(\*\*[^*]+\*\*|`[^`]+`)/g).map((part, i) => {
    if (part.startsWith("**") && part.endsWith("**") && part.length > 4) return <b key={i}><RichText text={part.slice(2, -2)} /></b>;
    if (part.startsWith("`") && part.endsWith("`") && part.length > 2) return <code key={i}>{part.slice(1, -1)}</code>;
    return <RichText key={i} text={part} />;
  });
}

export function Markdown({ text }: { text: string }) {
  const blocks: ReactNode[] = [];
  let list: string[] = [];
  let para: string[] = [];
  const flush = () => {
    if (para.length) { blocks.push(<p key={blocks.length}>{inline(para.join(" "))}</p>); para = []; }
    if (list.length) { blocks.push(<ul key={blocks.length}>{list.map((t, i) => <li key={i}>{inline(t)}</li>)}</ul>); list = []; }
  };
  for (const raw of text.split("\n")) {
    const line = raw.trimEnd();
    const head = /^(#{1,3})\s+(.*)$/.exec(line);
    if (head) {
      flush();
      const Tag = (`h${head[1].length}`) as "h1" | "h2" | "h3";
      blocks.push(<Tag key={blocks.length}>{inline(head[2])}</Tag>);
    } else if (/^\s*[-*]\s+/.test(line)) {
      if (para.length) flush();
      list.push(line.replace(/^\s*[-*]\s+/, ""));
    } else if (!line.trim()) flush();
    else {
      if (list.length) flush();
      para.push(line.trim());
    }
  }
  flush();
  return <div className="doc-body">{blocks}</div>;
}
