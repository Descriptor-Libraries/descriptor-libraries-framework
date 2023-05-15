import './App.css';
import { BrowserRouter as Router } from 'react-router-dom';

import ResponsiveAppBar from './common/Navbar';
import AppRoutes from './common/Routes';
import Footer from './common/Footer';

import Home from './pages/Home'
import About from './pages/About'
import Download from './pages/Download'
import Search  from './pages/SearchHook';
import NeighborSearch  from './pages/Neighbors';

const pages = {
  'Home': <Home />, 
  'About': <About />, 
  'Search': <Search />,
  'Neighbors': <NeighborSearch/>, 
  'Download': <Download />,
};

function App() {
  return (
    <div className="App">
      <Router>
        <header className="kraken-webapp">
          <ResponsiveAppBar pages={pages} />
        </header>
        <AppRoutes pages={pages} />
        <Footer />
      </Router>
    </div>
  );
}

export default App;
