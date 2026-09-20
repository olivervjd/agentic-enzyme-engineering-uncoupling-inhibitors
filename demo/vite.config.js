import { defineConfig } from 'vite';
import { fileURLToPath } from 'node:url';

export default defineConfig({
  server: { host: '127.0.0.1', fs: { strict: true, allow: [fileURLToPath(new URL('.', import.meta.url))] } },
  build: { rollupOptions: { output: { manualChunks: { three: ['three'], icons: ['lucide'] } } } }
});
