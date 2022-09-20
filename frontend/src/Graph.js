import React from 'react';
import Plot from 'react-plotly.js';

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
    fetch('/api/v1/molecule/molecules/umap?limit=500&category=pcn')
      .then(response => response.json())
      .then(items => this.setState({ pcn: items })
      );

    fetch('/api/v1/molecule/molecules/umap?limit=500&category=pc3')
      .then(response => response.json())
      .then(items => this.setState({ pc3: items })
      );

    fetch('/api/v1/molecule/molecules/umap?limit=500&category=po3')
      .then(response => response.json())
      .then(items => this.setState({ po3: items })
      );
  }

  render() {

    let myPlot = <Plot
    data={[
      {
        x: this.state.pcn.map( row => { return row.umap1 }),
        y: this.state.pcn.map( row => { return row.umap2 }),
        text: this.state.pcn.map( row => { return row.smiles }),
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'red'},
        name: 'PCN'
      },
      
      {
        x: this.state.pc3.map( row => { return row.umap1 }),
        y: this.state.pc3.map( row => { return row.umap2 }),
        //text: this.state.pc3.map( row => { return fetch(`/depict/bow/svg?smi=${row.smiles}&w=80&h=50&abbr=on&hdisp=bridgehead&showtitle=false&zoom=1.6&annotate=none`) } ),
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'blue'},
        name: 'PC3'
      },

      {
        x: this.state.po3.map( row => { return row.umap1 }),
        y: this.state.po3.map( row => { return row.umap2 }),
        text: this.state.po3.map( row => { return row.smiles }),
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'orange'},
        name: 'PO3'
      },
    ]}
    layout={ { 
      width: '100%', 
      height: '75%', 
      title: this.props.title,
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

  
    return (
      myPlot
    );
  }
}

export default MyPlot;
