// plotly.js-dist-min 没有自带类型；这里只按 any 用。
declare module "plotly.js-dist-min";

/** Vite asset URL imports (e.g. icons.svg?url). */
declare module "*?url" { const url: string; export default url; }

/** Vite build-time env (BASE_URL = the "base" in vite.config.ts). */
interface ImportMeta { readonly env: { readonly BASE_URL: string; readonly [key: string]: string | boolean | undefined } }

interface Window { mzdeskReady?: () => void }       // set by the splash script in index.html
