import React, { useEffect, useState } from 'react';
import Plot from 'react-plotly.js';
import { styled } from '@mui/material/styles';
import { TextField, Typography } from "@mui/material";
import Paper from '@mui/material/Paper';
import Grid from '@mui/material/Grid';
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import CircularProgress from '@mui/material/CircularProgress';

import Button from '@mui/material/Button';
import MenuItem from '@mui/material/MenuItem';
import FormGroup from '@mui/material/FormGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import Checkbox from '@mui/material/Checkbox';

const Item = styled(Paper)(({ theme }) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.text.secondary,
  }));

async function NeighborSearch(molecule_id, type, components, limit=48, skip=0) {
  /**
   * Requests neighbor data on the molecule from the backend.
   * @param {number} molecule_id Id of the molecule to query.
   * @param {string} type Type of dimensionality reduction. Can be one of PCA or UMAP.
   * @param {string} components String of comma separated integers.
   * @param {number} limit Limit of the search.
   * @param {number} skip How many values to skip.
   * @return {json}  The response json.
   */
    let encoded = encodeURIComponent(components);

    const response =  await fetch(`/api/molecules/${molecule_id}/neighbors/?type=${type}&components=${encoded}&skip=${skip}&limit=${limit}`)

    if (!response.ok) {
        throw new Error('Invalid Molecule Id')
    }

    else {
        return await response.json()
    }
}

async function retrieveSVG( smiles, distance ) {
  /**
   * Retrieves an svg for a molecule smiles.
   * @param {string} smiles Molecule smile representation.
   * @param {number} distance Distance between the molecule and the target.
   * @return {dictionary} A dictionary for each item which contains its svg, smiles string and distance from the target molecule.
   */
  let encoded = encodeURIComponent(smiles);

  const response = await fetch(`depict/cow/svg?smi=${encoded}&w=40&h=40`);
  
  const svg = await response.text();
  let result = {}
  result["svg"] = svg;
  result["smiles"] = smiles;
  result["distance"] = distance;
  return result
}

async function retrieveAllSVGs( items ) {
  /**
   * Retrieves all SVGs for the neighbor molecules.
   * @param {json} items Json containing all the molecule data returned from the neighbor search.
   * @return {json} Json containing the svg, smiles and distance for each neighbor.
   */
  return await Promise.all( items.map( (item) => { 
      return retrieveSVG(item.smiles, item.dist)
   } ) )
}

function dynamicGrid( svgs ) {
  /**
   * Generates a Grid filled with the svgs for the molecules.
   * @param {json} svgs Json containing the svg, smiles, and distance for each neighbor.
   * @return {jsx} The front end JSX that will generate the HTML users will interact with. It is a grid of items (described above) filled with the 
   * svg, smiles and distance of the molecule.
   */
  return (
      <Container>
      <Grid container spacing={2} sx= {{ mt: 3 }}>
      {
      svgs.map((result) => (
      <Grid item xs={12} md={4}>
          {result.distance == 0 ? <Item sx={{border: 3, borderColor: '#ed1c24'}}>
          <img alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(result.svg)}`} />
          <Typography sx={{ wordBreak: "break-word" }}> <strong>Smiles: </strong> { result.smiles }</Typography>
          </Item> :
          <Item>
          <img alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(result.svg)}`} />
          <Typography sx={{ wordBreak: "break-word" }}> <strong>Smiles: </strong> { result.smiles }</Typography>
          <Typography sx={{ wordBreak: "break-word" }}> <strong>Distance: </strong> { result.distance.toExponential(2) }</Typography>
          </Item>} 
      </Grid>
      ))
      }
      
  </Grid>
  </Container>
  )
}

export default function NeighborSearchHook () {
    /**
     * React functional component to generate the neighborhood UI.
     * @returns {jsx} The front end JSX that will generate the HTML users will interact with. It contains a search bar as well as generated 
     * a graph, and the dynamic grid component based on what data is available.
     */
    const interval = 15;
    const axis_dict = {"pca1": "pc1", "pca2": "pc2", "pca3": "pc3", "pca4": "pc4", "umap1": "umap1", "umap2": "umap2"};
    
    const [ moleculeid, setSearch ] = useState(1);
    const [ type, setType ] = useState("pca");
    const [ skip, setSkip ] = useState(0);
    const [ validMolecule, setValidMolecule ] = useState(true);
    const [ svg_results, setSVGResults ] = useState([])
    const [ searchPage, setSearchPage ] = useState(1);
    const [ isLoading, setIsLoading ] = useState(true);
    const [ searchToggle, setSearchToggle ] = useState(true);
    const [ isLoadingMore, setIsLoadingMore ] = useState(false);
    const [ molData, setMolData] = useState([]);
    const [ componentArrayForm, setComponentArrayForm ] = useState(["1", "2"]);

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

    function Graph(){
        /**
         * Generates plotly react graph. The target molecule is a red triangle facing up while the neighbors are grey triangles facing down.
         * The x and y labels are mapped to the axis dictionary to write pc1 instead of pca1.
         */
        // Shifting the data by 1, to avoid overwriting the target of the search
        let neighbors = molData.slice(1);
        let myPlot = <Plot onHover={ (event) => showSVG(event) } 
        onUnhover={ (event)=> hideSVG(event) } 
        style={{'width': '100%', 'height': '100%' }}
        useResizeHandler={true}
        data={[
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
            },
            // Creating the data series for the neighbors, by type
            // pc3
            {
              x: neighbors.map( row => {if (row.pat == "pc3") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "pc3") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
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
            // pn3
            {
              x: neighbors.map( row => {if (row.pat == "pn3") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "pn3") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
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
            // po3
            {
              x: neighbors.map( row => {if (row.pat == "po3") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "po3") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
              hovermode: "closest",
              type: 'scatter',
              mode: 'markers',
              marker: {color: 'yellow', size: 12 , 
                      symbol: "triangle-down",
                        line: {
                          width: 2,
                          color: 'DarkSlateGrey'}},
              name: 'P[O]<sub>3</sub>'
            },
            // pcn
            {
              x: neighbors.map( row => {if (row.pat == "pcn") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "pcn") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
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
            // phal
            {
              x: neighbors.map( row => {if (row.pat == "phal") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "phal") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
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
            // pon
            {
              x: neighbors.map( row => {if (row.pat == "pon") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "pon") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
              hovermode: "closest",
              type: 'scatter',
              mode: 'markers',
              marker: {color: 'orange', size: 12,
                      symbol: 'circle',
                      line: {
                        width: 2,
                        color: 'DarkSlateGrey'}},
              name: 'P[O]<sub>n</sub>[N]<sub>m</sub>'
            },
            // pco
            {
              x: neighbors.map( row => {if (row.pat == "pco") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "pco") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
              hovermode: "closest",
              type: 'scatter',
              mode: 'markers',
              marker: {color: 'purple', size: 12,
                      symbol: 'triangle-down',
                      line: {
                        width: 2,
                        color: 'DarkSlateGrey'}},
              name: 'P[C]<sub>n</sub>[O]<sub>m</sub>'
            },
            // psi
            {
              x: neighbors.map( row => {if (row.pat == "psi") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "psi") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
              hovermode: "closest",
              type: 'scatter',
              mode: 'markers',
              marker: {color: 'pink', size: 12,
                      symbol: 'hexagon2',
                      line: {
                        width: 2,
                        color: 'DarkSlateGrey'}},
              name: 'P[S]<sub>n</sub>[I]<sub>m</sub>'
            },
            // other
            {
              x: neighbors.map( row => {if (row.pat == "other") { return row.components[0] } }),
              y: neighbors.map( row => {if (row.pat == "other") { return row.components[1] } }),
              text: neighbors.map( row => { return encodeURIComponent(row.smiles) }),
              hovertemplate: "( %{x}, %{y})",
              hovermode: "closest",
              type: 'scatter',
              mode: 'markers',
              marker: {color: 'black', size: 12,
                      symbol: 'square',
                      line: {
                        width: 2,
                        color: 'DarkSlateGrey'}},
              name: 'Other'
            }
        ]}
        layout={ { 
          autosize: true,
          useResizeHandler: true,
          style: {width: '100%', height: '100%'},
          xaxis: {
            title: {
              text: axis_dict[molData[0].type + componentArrayForm[0]],
              font: {
                size: 18,
                color: '#7f7f7f'
              }
          }
        },
    
        yaxis: {
          title: {
            text: axis_dict[molData[0].type + componentArrayForm[1]],
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

    function buildComponentArray(event, label){
      /**
       * Manages the array of components based on how many checkboxes are ticked.
       * @param {event} event Selecting or deselecting the checkbox event.
       * @param {string} label The label of the checkbox.
       */

      // Add new label if the checkbox was checked
      if (event === true) {
        setComponentArrayForm(prevArray => {
          const newArray = [...prevArray, label];
          newArray.sort();
          return newArray;
        });
      }
      // Remove label otherwise
      else {
        setComponentArrayForm(prevArray => {
          const newArray = prevArray.filter(item => item !== label);
          newArray.sort();
          return newArray;
        });
      }
    }

    function arrayToString(components){
      /**
       * Converts components from an array to a comma separated string to use in our neighbor search API.
       * @param {string[]} components An array of components as strings.
       * @returns A comma separated string of the components.
       */
      return components.join(", ");
    }

    function loadMore() {
      /**
       * Loads more neighbors into the UI. 
       */
        setSkip(skip => skip + interval);
        setSearchPage( searchPage => searchPage + 1);
        setIsLoadingMore(true);
        loadNeighbors();
    }

    function newSearch() {
      /**
       * Searches for new neighbors. Resets alot of the props to their original state.
       */
        setSkip(0);
        setSearchPage(1);
        setSVGResults([]);
        // Just need to toggle this to make sure it toggles
        // so that effect will be triggered
        setIsLoading(true);
        setSearchToggle(!searchToggle);
        setMolData([]);
    }
 
    function loadNeighbors() {
      /**
       * Main driver function which loads the neighbors for a molecule requested by the user.
       * 
       */

        const fetchData = async () => {
            const molecule_data = await NeighborSearch(moleculeid, type, arrayToString(componentArrayForm), interval, skip);
            const svg_data = await retrieveAllSVGs(molecule_data);

            return [ molecule_data, svg_data ]
        }

        fetchData()
        .catch( (error) => {
            console.log(error) 
            setValidMolecule(false);
            setSVGResults([])
            setIsLoading(false)
            setIsLoadingMore(false) 
        } )
        .then( (items )=> {

            if (searchPage == 1) {
            setMolData(items[0]);
            setSVGResults(items[1]);
            }

            else {
                setMolData(molData.concat(items[0]));
                setSVGResults(svg_results.concat(items[1]));
            }

            setIsLoading(false);
            setIsLoadingMore(false);
            setValidMolecule(true);

          })

    }

    // initial load of data 
    useEffect( ( ) => { 
        loadNeighbors() }, 

        // eslint-disable-next-line react-hooks/exhaustive-deps
        [ searchPage, searchToggle ] 
    );
    
    // new search if dropdown changes
    useEffect(() => {
      newSearch();
      loadNeighbors();
    }, [type]);

    // new search if the molecule id changes
    useEffect(() => {
      newSearch();
      loadNeighbors();
    }, [moleculeid]);

    // new search if components change.
    useEffect(() => {
    newSearch();
    loadNeighbors();
  }, [componentArrayForm]);

    return (
        <Container maxWidth="lg">
        <h2>Neighbor Search</h2>
        <TextField 
                  style = {{width: 350}}
                  sx={{ m: 0.5}}
                  id="search-outline" 
                  label="Enter a Molecule ID to Search" 
                  variant="outlined"
                  defaultValue= {moleculeid} 
                  onChange = { event => setSearch( event.target.value ) }
        />
        <TextField
            sx={{ m: 0.5 }}
            select
            id="dimension-outline"
            value={type}
            // Need to clear MolData and SVGResults here. This is okay because when we switch type we need to do a new search right after.
            onChange={ function(event) {setType(event.target.value); setComponentArrayForm(["1", "2"]); setMolData([]); setSVGResults([]) }}
        >
            <MenuItem value={"pca"}>PCA</MenuItem>
            <MenuItem value={"umap"}>UMAP</MenuItem>
        </TextField>
        {type == "pca" && <FormGroup sx={{position: 'flex', flexDirection: 'row', justifyContent: 'center', alignItems:'center'}}>
                          <FormControlLabel control={<Checkbox defaultChecked value={"1"} onChange = {event => buildComponentArray(event.target.checked, event.target.value)}/>} label="1" />
                          <FormControlLabel control={<Checkbox defaultChecked value={"2"} onChange = {event => buildComponentArray(event.target.checked, event.target.value)}/>} label="2" />
                          <FormControlLabel control={<Checkbox onChange = {event => buildComponentArray(event.target.checked, "3")}/>} label="3" />
                          <FormControlLabel control={<Checkbox onChange = {event => buildComponentArray(event.target.checked, "4")}/>} label="4" />
                        </FormGroup>
        }
        {type == "umap" && <FormGroup sx={{position: 'flex', flexDirection: 'row', justifyContent: 'center', alignItems:'center'}}>
                          <FormControlLabel control={<Checkbox defaultChecked value={"1"} onChange = {event => buildComponentArray(event.target.checked, event.target.value)}/>} label="1" />
                          <FormControlLabel control={<Checkbox defaultChecked value={"2"} onChange = {event => buildComponentArray(event.target.checked, event.target.value)}/>} label="2" />
                        </FormGroup>
        }
        { isLoadingMore ? <CircularProgress sx={{ color: "#ed1c24" }} /> : <Button variant="contained" style={{backgroundColor: "#ed1c24"}} sx={{ m: 0.5 }} onClick={ () => loadMore() }>Load More</Button> }
        <Container sx={{justifyContent: 'center', my: 3}}>
            <Box sx={{ display: 'flex' }}>
            {/* If molecule is not valid and there is no mol data, then state that there are no results for the molecule ID requested*/}
            { !isLoading && !validMolecule && Object.keys(molData).length == 0 && <Typography>No results found for Molecule ID.</Typography> } 
            </Box>
            <Box>
            {/* If molecule is valid and there is mol data and the number of components is 2, then generate the graph based on the data*/}
            { !isLoading && validMolecule && Object.keys(molData).length > 0 && componentArrayForm.length == 2 && <Container sx={{ display: 'flex', height: 750}}>{ Graph() }</Container> } 
            </Box>
            <Box sx={{ display: 'flex' }}>
            {/* If molecule is valid and there is svg data, then generate the images of the molecules*/}
            { !isLoading && validMolecule && Object.keys(svg_results).length > 0 && 
             <Container> 
                { dynamicGrid(svg_results)  }
            </Container>  }
            
            </Box>
            { !isLoading && Object.keys(svg_results).length > 0 && 
            <Box sx={{ display: 'flex', justifyContent: 'center', mt: 2 }}>
            {isLoadingMore ? (
              <CircularProgress sx={{ color: "#ed1c24" }} />
            ) : (
              <Button
                variant="contained"
                style={{ backgroundColor: "#ed1c24" }}
                sx={{ m: 0.5 }}
                onClick={() => loadMore()}
              >
                Load More
              </Button>
            )}
          </Box>
          }

        </Container>
        </Container>
    )
}
