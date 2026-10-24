import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// 构建输出到 mzdesk/static，由 FastAPI 托管；开发时 /api 和 /lib 代理到后端。
// plotly.js 单独成块（约 1.5 MB gzip），第一次画图时才加载，所以放宽块大小警告。
export default defineConfig({
  plugins: [react()],
  build: { outDir: "../mzdesk/static", emptyOutDir: true, chunkSizeWarningLimit: 5000 },
  server: { port: 5173, proxy: { "/api": "http://localhost:8599", "/lib": "http://localhost:8599" } },
});
