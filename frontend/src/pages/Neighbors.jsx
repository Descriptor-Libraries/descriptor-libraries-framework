import React, { useEffect, useState } from 'react';
import { TextField, Typography } from "@mui/material";
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import CircularProgress from '@mui/material/CircularProgress';
import { ThemeProvider } from '@mui/material/styles';

import Button from '@mui/material/Button';
import MenuItem from '@mui/material/MenuItem';
import FormGroup from '@mui/material/FormGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import Checkbox from '@mui/material/Checkbox';

import Graph from '../components/Graph'

import { retrieveAllSVGs, dynamicGrid, theme } from '../common/MoleculeUtils';


async function NeighborSearch(molecule_id, type, components, limit=48, skip=0, signal) {
  /**
   * Requests neighbor data on the molecule from the backend.
   * @param {number} molecule_id Id of the molecule to query.
   * @param {string} type Type of dimensionality reduction. Can be one of PCA or UMAP.
   * @param {string} components String of comma separated integers.
   * @param {number} limit Limit of the search.
   * @param {number} skip How many values to skip.
   * @param {AbortSignal} signal Abortsignal object.
   * @return {json}  The response json.
   */
    let encoded = encodeURIComponent(components);

    const response =  await fetch(`/api/molecules/${molecule_id}/neighbors/?type=${type}&components=${encoded}&skip=${skip}&limit=${limit}`, {signal: signal})

    if (!response.ok) {
        throw new Error('Invalid Molecule Id')
    }

    else {
        return await response.json()
    }
}

export default function NeighborSearchHook () {
    /**
     * React functional component to generate the neighborhood UI.
     * @returns {jsx} The front end JSX that will generate the HTML users will interact with. It contains a search bar as well as generated 
     * a graph, and the dynamic grid component based on what data is available.
     */
    const interval = 15;
    
    const [ moleculeid, setSearch ] = useState(1);
    const [ type, setType ] = useState("pca");
    const [ graphType, setGraphType ] = useState("pca");
    const [ skip, setSkip ] = useState(0);
    const [ validMolecule, setValidMolecule ] = useState(true);
    const [ svg_results, setSVGResults ] = useState([])
    const [ searchPage, setSearchPage ] = useState(1);
    const [ isLoading, setIsLoading ] = useState(true);
    const [ isLoadingMore, setIsLoadingMore ] = useState(false);
    const [ molData, setMolData] = useState([]);
    const [ componentArrayForm, setComponentArrayForm ] = useState(["1", "2"]);
    const [ graphComponentArrayForm, setGraphComponentArrayForm ] = useState(["1", "2"]);
    const [ updatedParameters, setUpdatedParameters ] = useState(true);
    const [ searchToggle, setSearchToggle ] = useState(true);

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
        setSearchPage(searchPage => searchPage + 1);
        setIsLoadingMore(true);
        setSearchToggle(!searchToggle)
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
        setMolData([]);
        setSearchToggle(!searchToggle);
    }
 
    function loadNeighbors(signal) {
      /**
       * Main driver function which loads the neighbors for a molecule requested by the user.
       * @param {AbortSignal} signal Abortsignal object.
       */
        const fetchData = async () => {
            const molecule_data = await NeighborSearch(moleculeid, type, arrayToString(componentArrayForm), interval, skip, signal);
            const svg_data = await retrieveAllSVGs(molecule_data, signal);

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

    function updateGraphComponentProps() {
      /**
       * Used to update the type and component array form which is sent to the graph component. This is necessary because changing the type or components on the form
       * can change the output on the graph and its x and y axis immediately, even before a search is conducted. To prevent this, 2 sets of variables for type and components exist.
       * The type and componentArrayForm are used to hold the changes from the user while graphType and graphComponentArrayForm are used by the graph component to render the graph.
       * Thus this function is called to update the graphType and graphComponentArrayForm when a new search is called, so that the graph component can be updated, allowing the user to
       * make selections without updating the graph component.
       */
      setGraphType(type);
      setGraphComponentArrayForm(componentArrayForm);
    }

    // If any parameters change, we must set updatedParameters to true.
    useEffect(() => {
      setUpdatedParameters(true);
    }, [moleculeid, type, componentArrayForm])

    function _handleKeyDown(event) {
      if (event.key === "Enter") {
        newSearch();
        updateGraphComponentProps();
      }
    }

    // initial load of data
    // and load when search changes. 
    useEffect( ( ) => {
      const controller = new AbortController();
      const signal = controller.signal;

      setUpdatedParameters(false);
      loadNeighbors(signal);

      return () => {
        controller.abort();
      }
    },
      // eslint-disable-next-line react-hooks/exhaustive-deps
      [ searchToggle ]
    ); 

    return (
        <Container maxWidth="lg">
        <Typography variant="h2" textAlign="center">Neighbor Search</Typography>

        <Box display="flex" justifyContent="center" pt={2}>
        <TextField 
                  style = {{width: 350}}
                  sx={{ m: 0.5}}
                  id="search-outline" 
                  label="Enter a Molecule ID to Search" 
                  variant="outlined"
                  value= {moleculeid} 
                  onKeyDown = { (e) => _handleKeyDown(e) }
                  onChange = { event => setSearch( event.target.value ) }
                  InputProps={{endAdornment: <Button onClick={ () => {newSearch(); updateGraphComponentProps();} } >Search</Button>}}
        />
        <TextField
            sx={{ m: 0.5 }}
            select
            id="dimension-outline"
            value={type}
            onChange={ function(event) {setType(event.target.value); setComponentArrayForm(["1", "2"]);}}
        >
            <MenuItem value={"pca"}>PCA</MenuItem>
            <MenuItem value={"umap"}>UMAP</MenuItem>
        </TextField>
        </Box>
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
        <Box display="flex" justifyContent="center">
        { (isLoading || isLoadingMore) ? <CircularProgress />: <ThemeProvider theme={theme}> <Button disabled={updatedParameters} variant="contained" sx={{ m: 0.5 }} onClick={ () => loadMore() }>Load More</Button> </ThemeProvider>}
        </Box>
        <Container sx={{justifyContent: 'center', my: 3}}>
            <Box sx={{ display: 'flex' }} justifyContent="center">
            {/* If molecule is not valid and there is no mol data, then state that there are no results for the molecule ID requested*/}
            { !isLoading && !validMolecule && Object.keys(molData).length == 0 && <Typography>No results found for Molecule ID.</Typography> } 
            </Box>
            <Box>
            {/* If molecule is valid and there is mol data, then generate the graph based on the data*/}
            { !isLoading && validMolecule && Object.keys(molData).length > 0 && graphComponentArrayForm.length > 1 && <Container sx={{ display: 'flex', height: 750}}>{ <Graph molData={molData} componentArray={graphComponentArrayForm} type={graphType} neighborSearch={true}></Graph> }</Container> } 
            </Box>
            <Box sx={{ display: 'flex' }} justifyContent="center">
            {/* If molecule is valid and there is svg data, then generate the images of the molecules*/}
            { !isLoading && validMolecule && Object.keys(svg_results).length > 0 && 
             <Container> 
                { dynamicGrid(svg_results)  }
            </Container>  }
            
            </Box>
            { !isLoading && Object.keys(svg_results).length > 0 && 
            <Box sx={{ display: 'flex', justifyContent: 'center', mt: 2 }}>
            {isLoadingMore ? (
              <CircularProgress />
            ) : (
              <ThemeProvider theme={theme}> <Button disabled={updatedParameters} variant="contained" sx={{ m: 0.5 }} onClick={ () => loadMore() }>Load More</Button> </ThemeProvider>
            )}
          </Box>
          }

        </Container>
        </Container>
    )
}
