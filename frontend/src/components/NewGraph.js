import React from 'react';
import Plot from 'react-plotly.js';

export default function Graph({ molData, componentArray, neighborSearch }){

    const axis_dict = {"pca1": "pc1", "pca2": "pc2", "pca3": "pc3", "pca4": "pc4", "umap1": "umap1", "umap2": "umap2"};
    
    // Plotting functions to show molecules on hover
    function showSVGWindow(svg, event) {
    /**
     * Creates SVG window on the molecule element when hovering of a point on the plotly graph.
     * @param {string} svg Svg of the molecule.
     * @param {event} event Hover event when hovering over a point on the plotly graph.
     */

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
        mol.style.position = "fixed";
        mol.style.left = `${xpos}px`;
        mol.style.top = `${ypos}px`;
        }
    
    function showSVG(event) {
    /**
     * Requests svg data for the molecule you are hovering on.
     * @param {event} event Hover even when hovering over a point on the plotly graph.
     */
    fetch(`depict/cow/svg?smi=${event.points[0].text}&w=40&h=40`).then(response => 
        response.text() ).then( body => showSVGWindow(body, event) );
    }
    
    function hideSVG() {
    /**
     * Removes molecule element which holds its SVG, when you are no longer hovering.
     */
    if (document.getElementById("molecule")) {
        document.getElementById("molecule").remove()
    }
    }
    
    /**
     * Generates plotly react graph.
     * The x and y labels are mapped to the axis dictionary to write pc1 instead of pca1.
     */
    // Shifting the data by 1, to avoid overwriting the target of the search if we need it.
    let values;
    if (neighborSearch) {
        values = molData.slice(1);
    }
    else {
        values = molData;
    }

    let myPlot = <Plot 
                    onHover={ (event) => showSVG(event) } 
                    onUnhover={ (event)=> hideSVG(event) } 
                    style={{'width': '100%', 'height': '100%' }}
                    useResizeHandler={true}

                    // Starting data generation
                    data={[
                        // Creating the data series for the values, by type
                        // pc3
                        {
                            x: values.map( row => {if (row.pat == "pc3") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "pc3") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'SlateGrey', size: 12 ,
                                    symbol: 'triangle-down', 
                                    line: {
                                        width: 2,
                                        color: 'DarkSlateGrey'}},
                            name: 'P[C]<sub>3</sub>',
                            showlegend: values.some((element) => element.pat == "pc3")
                        },
                        // pn3
                        {
                            x: values.map( row => {if (row.pat == "pn3") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "pn3") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'blue', size: 12,
                                    symbol: 'circle',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'P[N]<sub>3</sub>',
                            showlegend: values.some((element) => element.pat == "pn3")
                        },
                        // po3
                        {
                            x: values.map( row => {if (row.pat == "po3") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "po3") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'yellow', size: 12 , 
                                    symbol: "triangle-down",
                                    line: {
                                        width: 2,
                                        color: 'DarkSlateGrey'}},
                            name: 'P[O]<sub>3</sub>',
                            showlegend: values.some((element) => element.pat == "po3")
                        },
                        // pcn
                        {
                            x: values.map( row => {if (row.pat == "pcn") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "pcn") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'DarkBlue', size: 12,
                                    symbol: 'hexagon',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'P[C]<sub>n</sub>[N]<sub>m</sub>',
                            showlegend: values.some((element) => element.pat == "pcn")
                        },
                        // phal
                        {
                            x: values.map( row => {if (row.pat == "phal") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "phal") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'LimeGreen', size: 12,
                                    symbol: 'hexagon2',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'PF<sub>n</sub>[R]<sub>m</sub>',
                            showlegend: values.some((element) => element.pat == "phal")
                        },
                        // pon
                        {
                            x: values.map( row => {if (row.pat == "pon") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "pon") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'orange', size: 12,
                                    symbol: 'circle',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'P[O]<sub>n</sub>[N]<sub>m</sub>',
                            showlegend: values.some((element) => element.pat == "pon")
                        },
                        // pco
                        {
                            x: values.map( row => {if (row.pat == "pco") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "pco") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'purple', size: 12,
                                    symbol: 'triangle-down',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'P[C]<sub>n</sub>[O]<sub>m</sub>',
                            showlegend: values.some((element) => element.pat == "pco")
                        },
                        // psi
                        {
                            x: values.map( row => {if (row.pat == "psi") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "psi") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'pink', size: 12,
                                    symbol: 'hexagon2',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'P[S]<sub>n</sub>[I]<sub>m</sub>',
                            showlegend: values.some((element) => element.pat == "psi")
                        },
                        // other
                        {
                            x: values.map( row => {if (row.pat == "other") { return row.components[0] } }),
                            y: values.map( row => {if (row.pat == "other") { return row.components[1] } }),
                            text: values.map( row => { return encodeURIComponent(row.smiles) }),
                            hovertemplate: "( %{x}, %{y})",
                            hovermode: "closest",
                            type: 'scatter',
                            mode: 'markers',
                            marker: {color: 'black', size: 12,
                                    symbol: 'square',
                                    line: {
                                    width: 2,
                                    color: 'DarkSlateGrey'}},
                            name: 'Other',
                            showlegend: values.some((element) => element.pat == "other")
                        }
                    ]}
                    layout={ { 
                        autosize: true,
                        useResizeHandler: true,
                        style: {width: '100%', height: '100%'},
                        xaxis: {
                        title: {
                            text: axis_dict[molData[0].type + componentArray[0]],
                            font: {
                            size: 18,
                            color: '#7f7f7f'
                            }
                        }
                    },

                    yaxis: {
                        title: {
                        text: axis_dict[molData[0].type + componentArray[1]],
                        font: {
                            size: 18,
                            color: '#7f7f7f'
                        }
                    }
                    }
                    } }
                />
    
    if (neighborSearch) {
        myPlot.props.data.push(                        
        // Creating the data series for the target of the search
        {
        x: [molData[0].components[0]],
        y: [molData[0].components[1]],
        text: [encodeURIComponent(molData[0].smiles)],
        hovertemplate: "( %{x}, %{y})",
        hovermode: "closest",
        type: 'scatter',
        mode: 'markers',
        marker: {color: 'red', size: 12 , 
                symbol: "triangle-up",
                line: {
                    width: 2,
                    color: 'DarkSlateGrey'}},
        name: 'Target'
        },)
    }

    return (
        myPlot
        );
}