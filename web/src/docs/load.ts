import { BASE } from "../api";

/** The text of a Docs page (GET /api/docs/<id>). null = there is no file for it (404); any other failure throws. */
export async function loadPage(id: string): Promise<string | null> {
  const r = await fetch(`${BASE}api/docs/${encodeURIComponent(id)}`);
  if (r.status === 404) return null;
  if (!r.ok) throw new Error(`${r.status} ${r.statusText}`);
  return (await r.json()).markdown as string;
}
