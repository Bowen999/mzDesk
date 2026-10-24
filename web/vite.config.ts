import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

// 构建输出到 mzdesk/static，由 FastAPI 托管；开发时 /api 代理到后端。
export default defineConfig({
  plugins: [react()],
  build: { outDir: "../mzdesk/static", emptyOutDir: true },
  server: { port: 5173, proxy: { "/api": "http://localhost:8599" } },
});
