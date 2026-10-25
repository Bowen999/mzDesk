// plotly.js-dist-min 没有自带类型；这里只按 any 用。
declare module "plotly.js-dist-min";

/** Vite asset URL imports (e.g. icons.svg?url). */
declare module "*?url" { const url: string; export default url; }
