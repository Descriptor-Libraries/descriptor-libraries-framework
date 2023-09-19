import React, { lazy, Suspense } from 'react';
import { BrowserRouter as Router, useLocation } from 'react-router-dom';

import ResponsiveAppBar from './common/Navbar';
import { AppRoutes } from './common/Routes';
import Footer from './common/Footer';
import { Box } from '@mui/material'

import Home from './pages/Home';

const Search = lazy(() => import('./pages/SearchHook'));
const NeighborSearch = lazy(() => import('./pages/Neighbors'));
const About = lazy(() => import('./pages/About'))

const pages = {
  'Home': <Home />, 
  'Search': <Search />,
  'Neighbors': <NeighborSearch/>, 
  'About': <About/>,
};

function Content() {
  const location = useLocation();
  return (
    <>
      {location.pathname !== "/" && location.pathname !== "/home" && <ResponsiveAppBar pages={pages} />}
      <Suspense fallback={<div>Loading...</div>}>
        <AppRoutes pages={pages} />
      </Suspense>
    </>
  );
}

function App() {
  return (
    <>
    <Box sx={{ display: 'flex', flexDirection: 'column'}}>
        <Router>
            <Content />
        </Router>
    </Box>
    <Footer />
    </>
  );
}

export default App;
