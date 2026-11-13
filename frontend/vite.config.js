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
      // Ketcher reads `process.env.*` in the browser. Give it only the values it
      // needs. Never pass the whole `process.env` here: that would copy every
      // environment variable of the build machine into the public bundle.
      'process.env': JSON.stringify({
        NODE_ENV: mode === 'production' ? 'production' : 'development',
        PUBLIC_URL: '',
      }),
    },
    build: {
      commonjsOptions: {
        transformMixedEsModules: true,
      },
    },
  }
});