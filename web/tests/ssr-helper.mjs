// Shared setup of the component tests that render with react-dom/server: a Vite SSR server over src/, a render helper, and a quiet console
// (the "?" popovers use useLayoutEffect, which React warns about on the server only; nothing there is wrong).
import { mkdtempSync, rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import React from "react";
import { renderToStaticMarkup } from "react-dom/server";
import { createServer } from "vite";

export { React };

export async function startSsr(label) {
  // Figure's theme module initializes the document theme on import.
  globalThis.document = { documentElement: { removeAttribute() {}, setAttribute() {} } };
  const root = dirname(dirname(fileURLToPath(import.meta.url)));
  const cache = mkdtempSync(join(tmpdir(), `mzdesk-${label}-test-`));
  const server = await createServer({ configFile: false, root, cacheDir: cache, server: { middlewareMode: true, hmr: false }, appType: "custom" });
  const consoleError = console.error;
  console.error = (...args) => { if (!String(args[0]).includes("useLayoutEffect does nothing on the server")) consoleError(...args); };
  return {
    load: path => server.ssrLoadModule(path),
    render: (component, props) => renderToStaticMarkup(React.createElement(component, props)),
    close: async () => { console.error = consoleError; await server.close(); rmSync(cache, { recursive: true, force: true }); delete globalThis.document; },
  };
}
