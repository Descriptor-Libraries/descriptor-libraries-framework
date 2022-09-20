import logo from './logo.svg';
import './App.css';

import MyPlot from "./Graph.js";

function SmilesDepict(smiles) {
  console.log("hello")
}

function App() {
  return (
    <div className="App">
      <header className="App-header">
        umap plot
        <MyPlot 
          onHover={SmilesDepict}
          title="umap plot!"
          ></MyPlot>
      </header>
    </div>
  );
}

export default App;
