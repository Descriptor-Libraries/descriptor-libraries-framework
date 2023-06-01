import './App.css';
import { BrowserRouter as Router, useLocation} from 'react-router-dom';

import ResponsiveAppBar from './common/Navbar';
import { AppRoutes } from './common/Routes';
import Footer from './common/Footer';

import Home from './pages/Home'
import About from './pages/About'
import Download from './pages/Download'
import Search  from './pages/SearchHook';
import NeighborSearch  from './pages/Neighbors';

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
      <AppRoutes pages={pages} />
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
