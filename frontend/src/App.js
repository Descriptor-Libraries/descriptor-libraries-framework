import logo from './logo.svg';
import './App.css';

import MyPlot from "./Graph.js";
import ResponsiveAppBar from './Navbar';

function App() {
  return (
    <div className="App">
      <header className="kraken-webapp">
      </header>
      <ResponsiveAppBar />
      <MyPlot />
    </div>
  );
}

export default App;
