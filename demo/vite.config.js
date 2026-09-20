import { defineConfig } from 'vite';
import { fileURLToPath } from 'node:url';
export default defineConfig({
  base:'./',
  server:{host:'127.0.0.1',fs:{strict:true,allow:[fileURLToPath(new URL('.',import.meta.url))]}},
  build:{rollupOptions:{output:{inlineDynamicImports:true}}}
});
