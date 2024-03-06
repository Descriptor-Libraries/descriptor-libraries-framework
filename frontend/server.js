const express = require('express');
const path = require('path');
const app = express();

const PORT = process.env.PORT || 3000;
// Assuming VITE_BASE_URL is something like '/myapp/'
// Make sure it falls back to '/' if VITE_BASE_URL is not set.
const BASE_URL = process.env.VITE_BASE_URL || '/';

// Serve files from the Vite build output directory with the base URL

// Mounted public directory - comes first so it takes precedence  
app.use(BASE_URL, express.static(path.join(__dirname, 'public')));

// Pre-built dist directory.
app.use(BASE_URL, express.static(path.join(__dirname, 'dist')));

// Handle any requests that don't match the ones above by always returning the main index.html file
app.get('*', (req, res) => {
  res.sendFile(path.join(__dirname, 'dist', 'index.html'));
});

app.listen(PORT, () => {
    console.log(`Server running on http://localhost:${PORT}${BASE_URL}`);
});
