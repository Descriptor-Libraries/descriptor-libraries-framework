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
      pn3: [],
      phal: [],
    };
  }

  componentDidMount() {
    fetch('/api/molecules/umap?category=pcn')
      .then(response => response.json())
      .then(items => this.setState({ 
        pcn: items
      })
      );

    fetch('/api/molecules/umap?category=pc3')
      .then(response => response.json())
      .then(items => this.setState({ pc3: items })
      );

    fetch('/api/molecules/umap?category=po3')
      .then(response => response.json())
      .then(items => this.setState({ po3: items })
      );

      fetch('/api/molecules/umap?category=phal')
      .then(response => response.json())
      .then(items => this.setState({ phal: items })
      );

      fetch('/api/molecules/umap?category=pn3')
      .then(response => response.json())
      .then(items => this.setState({ pn3: items })
      );

  }

  render() {

    
    let myPlot = <div id="plotly-container" style={{'width': '60%', 'height': '85%', 'margin': 'auto' }}> <Plot onHover={ (event) => showSVG(event) } 
    onUnhover={ (event)=> hideSVG(event) } 
    style={{'width': '100%', 'height': '100%' }}
    useResizeHandler={true}
    data={[

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
        name: 'P[C]<sub>3</sub>'
      },

      {
        x: this.state.pn3.map( row => { return row.umap1 }),
        y: this.state.pn3.map( row => { return row.umap2 }),
        text: this.state.pn3.map( row => { return encodeURIComponent(row.smiles) }),
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
        x: this.state.po3.map( row => { return row.umap1 }),
        y: this.state.po3.map( row => { return row.umap2 }),
        text: this.state.po3.map( row => { return encodeURIComponent(row.smiles) }),
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
        x: this.state.pcn.map( row => { return row.umap1 }),
        y: this.state.pcn.map( row => { return row.umap2 }),
        text: this.state.pcn.map( row => { return encodeURIComponent(row.smiles) }),
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
        x: this.state.phal.map( row => { return row.umap1 }),
        y: this.state.phal.map( row => { return row.umap2 }),
        text: this.state.phal.map( row => { return encodeURIComponent(row.smiles) }),
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
