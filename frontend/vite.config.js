import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import { nodePolyfills } from 'vite-plugin-node-polyfills'; // have to do this for ketcher
import vitePluginRaw from 'vite-plugin-raw';

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  return {
    base: process.env.VITE_BASE_URL || '/base_url', 
    plugins: [
      react(), 
      nodePolyfills({
        include: ['buffer', 'process'],
        globals: {
          Buffer: true,
          global: true,
          process: true,
        },
      }),
      vitePluginRaw({
        match: /\.sdf|\.ket/,
      }),
    ],
    server: {
      host: '0.0.0.0',
      port: 3000,
    },
    define: {
      'process.env': process.env // have to do this for ketcher 
    },
    build: {
      commonjsOptions: {
        transformMixedEsModules: true,
      },
    },
  }
});