import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// ketcher-core 的 raphael-ext 在 ESM 里写了 require("raphael")，浏览器里没有 require：改成 import
const ketcherRaphael = {
  name: "ketcher-raphael-import",
  enforce: "pre" as const,
  transform(code: string, id: string) {
    if (!id.includes("raphael-ext.modern.js")) return null;
    return `import __raphael from "raphael";\n` + code.replace(/require\((['"])raphael\1\)/, "__raphael");
  },
};

// 构建输出到 mzdesk/static，由 FastAPI 托管；开发时 /api 和 /lib 代理到后端。
// plotly.js 单独成块（约 1.5 MB gzip），第一次画图时才加载，所以放宽块大小警告。
export default defineConfig({
  plugins: [ketcherRaphael, react()],
  // Ketcher 引用 Node 的 events 模块；浏览器里用同名 npm 包顶替
  resolve: { alias: { events: "events/" } },
  define: { "process.env": {}, global: "globalThis" },
  build: { outDir: "../mzdesk/static", emptyOutDir: true, chunkSizeWarningLimit: 5000 },
  server: { port: 5173, proxy: { "/api": "http://localhost:8599", "/lib": "http://localhost:8599" } },
});
