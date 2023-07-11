import React, { lazy, Suspense } from 'react';
import { BrowserRouter as Router, useLocation } from 'react-router-dom';

import ResponsiveAppBar from './common/Navbar';
import { AppRoutes } from './common/Routes';
import Footer from './common/Footer';

import Home from './pages/Home'
//import Download from './pages/Download'
//import Search from './pages/SearchHook'
//import NeighborSearch from './pages/Neighbors'

// Lazy loading for improved performance
const Download = lazy(() => import('./pages/Download'));
const Search = lazy(() => import('./pages/SearchHook'));
const NeighborSearch = lazy(() => import('./pages/Neighbors'));

const pages = {
  'Home': <Home />, 
  'Search': <Search />,
  'Neighbors': <NeighborSearch/>, 
  'Download': <Download />,
};

function Content() {
  const location = useLocation();
  return (
    <div class="content">
      {location.pathname !== "/" && location.pathname !== "/home" && <ResponsiveAppBar pages={pages} />}
      <Suspense fallback={<div>Loading...</div>}>
        <AppRoutes pages={pages} />
      </Suspense>
    </div>
  );
}


function App() {
  return (
    <div className="App">
      <Router>
        <header className="kraken-webapp">
          <Content />
        </header>
        <Footer />
      </Router>
    </div>
  );
}


export default App;
