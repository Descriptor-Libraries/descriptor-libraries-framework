import React from 'react';
import Plot from 'react-plotly.js';


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

    this.state = {
      pcn: [],
      pc3: [],
      po3: [],
    };
  }

  componentDidMount() {
    fetch('/api/v1/molecule/molecules/umap?category=pcn')
      .then(response => response.json())
      .then(items => this.setState({ 
        pcn: items
      })
      );

    fetch('/api/v1/molecule/molecules/umap?category=pc3')
      .then(response => response.json())
      .then(items => this.setState({ pc3: items })
      );

    fetch('/api/v1/molecule/molecules/umap?category=po3')
      .then(response => response.json())
      .then(items => this.setState({ po3: items })
      );
  }

  render() {

    
    let myPlot = <div id="plotly-container" style={{'width': '60%', 'height': '85%', 'margin': 'auto' }}> <Plot onHover={ (event) => showSVG(event) } 
    onUnhover={ (event)=> hideSVG(event) } 
    style={{'width': '100%', 'height': '100%' }}
    useResizeHandler={true}
    data={[
      {
        x: this.state.pcn.map( row => { return row.umap1 }),
        y: this.state.pcn.map( row => { return row.umap2 }),
        text: this.state.pcn.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'blue', size: 12,
          symbol: 'hexagon',
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'PCN'
      },
      
      {
        x: this.state.pc3.map( row => { return row.umap1 }),
        y: this.state.pc3.map( row => { return row.umap2 }),
        text: this.state.pc3.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'SlateGrey', size: 12 ,
        symbol: 'triangle-down', 
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'PC<sub>3</sub>'
      },

      {
        x: this.state.po3.map( row => { return row.umap1 }),
        y: this.state.po3.map( row => { return row.umap2 }),
        text: this.state.po3.map( row => { return encodeURIComponent(row.smiles) }),
        hovertemplate: "( %{x}, %{y} )",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'red', size: 12 , 
          line: {
            width: 2,
            color: 'DarkSlateGrey'}},
        name: 'PO<sub>3</sub>'
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
