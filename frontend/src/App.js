import logo from './logo.svg';
import './App.css';

import MyPlot from "./Graph.js";

function App() {
  return (
    <div className="App">
      <header className="App-header">
        umap plot
        <MyPlot title="umap plot!"
          ></MyPlot>
      </header>
    </div>
  );
}

export default App;
