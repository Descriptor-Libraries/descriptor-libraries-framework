import React from "react";
//import ReactDOM from "react-dom";
import { createRoot } from "react-dom/client"
import "./index.css";
import App from "./App";
import reportWebVitals from "./reportWebVitals";

function setLinkIcon() {
  const link = document.createElement('link');
  link.rel = 'icon';
  link.href = `${process.env.VITE_BASE_URL}/brand/favicon.ico`;
  document.head.appendChild(link);

}

function setSiteManifest() {
  const link = document.createElement('link');
  link.rel = 'manifest';
  link.href = `${process.env.VITE_BASE_URL}/manifest.json`;
}

function setAppBrand() {
  fetch(`${process.env.VITE_BASE_URL}/brand/names.json`)
  .then(response => response.json()) 
  .then(data => {
    // Update the site title
    document.title = data[0].title;

    // Update the site favicon
    setLinkIcon();

    // Update the site manifest
    setSiteManifest();
  })
  .catch(error => console.error('Error loading site configuration:', error));
}


const container = document.getElementById('root');

const root = createRoot(container)

setAppBrand();

root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>,
);

// If you want to start measuring performance in your app, pass a function
// to log results (for example: reportWebVitals(console.log))
// or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
reportWebVitals();