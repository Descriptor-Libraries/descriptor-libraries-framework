import React from 'react';
import Plot from 'react-plotly.js';

const datasetsMap = new Map([
  ['showPCN', '/api/molecules/umap?category=pcn'],
  ['showPC3', '/api/molecules/umap?category=pc3'],
  ['showPO3', '/api/molecules/umap?category=po3'],
  ['showPN3', '/api/molecules/umap?category=pn3'],
  ['showPHAL', '/api/molecules/umap?category=phal'],
]);


function showSVGWindow(svg, event) {

  // remove in case existing
  if (document.getElementById("molecule")) {
    document.getElementById("molecule").remove()
  }

  let mol = document.createElement("g")
  mol.setAttribute("id", "molecule")
  
  let plotly_container = document.getElementsByClassName("plotly")
  mol.innerHTML = svg;
  
  let xpos = event.event.clientX - 160;
  let ypos = event.event.clientY - 160;

  plotly_container[0].appendChild(mol);
  mol.style.position = "absolute";
  mol.style.left = `${xpos}px`;
  mol.style.top = `${ypos}px`;

}

function showSVG(event) {
  fetch(`depict/cow/svg?smi=${event.points[0].text}&w=40&h=40`).then(response => 
    response.text() ).then( body => showSVGWindow(body, event) );
}

function hideSVG(event) {
  if (document.getElementById("molecule")) {
    document.getElementById("molecule").remove()
  }

}

class MyPlot extends React.Component {

  // Constructor
  constructor(props) {
    super(props);
  
    const datasets = {};
    for (const [key, url] of datasetsMap) {
      const datasetKey = key.replace(/^show/, '').toLowerCase();
      datasets[datasetKey] = [];
    }
  
    this.state = { datasets };
  }
  

  componentDidMount() {
    console.log("load data")

    this.loadData()
     
  }

  componentDidUpdate(prevProps) {
    const { checkboxState } = this.props;

    if (Object.keys(checkboxState).some((key) => checkboxState[key] !== prevProps.checkboxState[key])) {
      this.loadData();
    }
  }

  loadData() {
    const { checkboxState } = this.props;
    const datasetsToLoad = [];

    for (const [key, url] of datasetsMap) {
      if (checkboxState[key]) {
        datasetsToLoad.push(fetch(url).then((response) => response.json()));
      }
    }

    console.log(datasetsToLoad)

    Promise.all(datasetsToLoad)
      .then(datasets => {
        this.setState({
          datasets: {
            pcn: checkboxState.showPCN ? datasets.shift() : [],
            pc3: checkboxState.showPC3 ? datasets.shift() : [],
            po3: checkboxState.showPO3 ? datasets.shift() : [],
            pn3: checkboxState.showPN3 ? datasets.shift() : [],
            phal: checkboxState.showPHAL ? datasets.shift() : [],
          }
        });
        console.log('datasets')
        console.log(this.datasets)
      })
      .catch(error => console.error('Error loading data:', error));
  }


  

  printMessage() {
    if (this.props.showPC3) {
      console.log('showPC3 is true');
    } else {
      console.log('showPC3 is false');
    }
  }

  render() {

    console.log('render is running')
    let myPlot = <div id="plotly-container" style={{'width': '60%', 'height': '85%', 'margin': 'auto' }}> <Plot onHover={ (event) => showSVG(event) } 
    onUnhover={ (event)=> hideSVG(event) } 
    style={{'width': '100%', 'height': '100%' }}
    useResizeHandler={true}
    data={[
       {
        x: this.state.datasets.pc3.map( row => { return row.umap1 }),
        y: this.state.datasets.pc3.map( row => { return row.umap2 }),
        text: this.state.datasets.pc3.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'SlateGrey', size: 12 ,
        symbol: 'triangle-down', 
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'P[C]<sub>3</sub>'
      },

      {
        x: this.state.datasets.pn3.map( row => { return row.umap1 }),
        y: this.state.datasets.pn3.map( row => { return row.umap2 }),
        text: this.state.datasets.pn3.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'blue', size: 12,
          symbol: 'circle',
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'P[N]<sub>3</sub>'
      },


      {
        x: this.state.datasets.po3.map( row => { return row.umap1 }),
        y: this.state.datasets.po3.map( row => { return row.umap2 }),
        text: this.state.datasets.po3.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'red', size: 12 , 
        symbol: "triangle-up",
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'P[O]<sub>3</sub>'
      },


      {
        x: this.state.datasets.pcn.map( row => { return row.umap1 }),
        y: this.state.datasets.pcn.map( row => { return row.umap2 }),
        text: this.state.datasets.pcn.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'DarkBlue', size: 12,
          symbol: 'hexagon',
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'P[C]<sub>n</sub>[N]<sub>m</sub>'
      },

      
      {
        x: this.state.datasets.phal.map( row => { return row.umap1 }),
        y: this.state.datasets.phal.map( row => { return row.umap2 }),
        text: this.state.datasets.phal.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'LimeGreen', size: 12,
          symbol: 'hexagon2',
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'PF<sub>n</sub>[R]<sub>m</sub>'
      },
      
    ]}
    layout={ { 
      title: this.props.title,
      autosize: true,
      useResizeHandler: true,
      style: {width: '100%', height: '100%'},
      xaxis: {
        title: {
          text: 'umap1',
          font: {
            size: 18,
            color: '#7f7f7f'
          }
      }
    },

    yaxis: {
      title: {
        text: 'umap2',
        font: {
          size: 18,
          color: '#7f7f7f'
        }
    }
  }
    
    } }
  />
</div>
  
    return (
      myPlot
    );
  }
}

export default MyPlot;
