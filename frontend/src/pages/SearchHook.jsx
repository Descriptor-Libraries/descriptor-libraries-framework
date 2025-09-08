import React, { useEffect, useState, useCallback } from 'react';
import { TextField, Typography } from "@mui/material";
import Grid from '@mui/material/Grid';
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import CircularProgress from '@mui/material/CircularProgress';

import { Switch } from '@mui/material';
import Button from '@mui/material/Button';

import FullScreenDialog from '../components/KetcherPopup';

import { substructureSearch, retrieveAllSVGs, dynamicGrid, extractIdsFromResults, downloadMoleculeData } from '../common/MoleculeUtils';

import DropDownButton from '../components/DataDownloadButton';

export default function SearchHook () {

    const interval = 15;

    const [ searchString, setSearch ] = useState('');
    const [ skip, setSkip ] = useState(0);
    const [ results, setResults ] = useState([]);
    const [ validSmiles, setValidSmiles ] = useState(true);
    const [ svg_results, setSVGResults ] = useState([])
    const [ searchPage, setSearchPage ] = useState(1);
    const [ isLoading, setIsLoading ] = useState(false);
    const [ searchToggle, setSearchToggle ] = useState(true);
    const [ isLoadingMore, setIsLoadingMore ] = useState(false);
    const [ smiles, setSmiles ] = useState('PC=C');
    const [ SMARTS, setSMARTS ] = useState('[#15]-[#6]=[#6]');
    const [ representation, setRepresentation ] = useState("smiles");
    const [ toggleRepresentation, setToggleRepresentation ] = useState(true);
    const [ switchCheck, setSwitchCheck ] = useState(true);
    const [ fromKetcher, setFromKetcher ] = useState(false);
    const [ updatedParameters, setUpdatedParameters ] = useState(false);
    const [ moleculeIDs, setMoleculeIDs ] = useState("");
    const [ availableDataTypes, setAvailableDataTypes ] = useState([]);

    const reverseMapping = {
      "ml_data": "ML Data",
      "dft_data": "DFT Data",
      "xtb_data": "xTB Data",
      "xtb_ni_data": "xTB_Ni Data"
  };
  
  
    useEffect(() => {
      fetch(`/api/${document.location.pathname.split('/')[1]}/molecules/data_types`)
      .then(response => response.json())
      .then(data => {
        const translatedDataTypes = data["available_types"].map(key => reverseMapping[key] || key);
        setAvailableDataTypes(translatedDataTypes);
      });
  }, []);

  // Load default SMILES from branding configuration
  useEffect(() => {
    fetch(`/${document.location.pathname.split('/')[1]}/brand/names.json`)
      .then(response => response.json())
      .then(data => {
        if (data && data[0] && data[0].default_smiles) {
          const defaultSmiles = data[0].default_smiles;
          setSearch(defaultSmiles);
          setSmiles(defaultSmiles);
        } else {
          // Fallback to hardcoded default if branding fails
          setSearch('C=C');
          setSmiles('C=C');
        }
      })
      .catch(error => {
        console.log('Could not load default SMILES from branding config:', error);
        // Fallback to hardcoded default
        setSearch('C=C');
        setSmiles('C=C');
      });
  }, []);

    // Extract the molecule ids from the results
    useEffect(() => {
      setMoleculeIDs(extractIdsFromResults(svg_results));
    }, [svg_results]);
  

    // Call back function to get the smiles and SMARTS from ketcher
    const ketcherCallBack = (newState) => {
        console.log('Ketcher callback received:', newState);
        console.log('Current representation:', representation);
        console.log('Current searchString before update:', searchString);

        if (newState[0] && newState[0].trim() !== '') {
          // Set the smiles and SMARTS for the current molecule
          setSmiles(newState[0]);
          setSMARTS(newState[1]);

          // This came from ketcher
          setFromKetcher(true);
          
          if (representation === "smiles"){
              console.log('Setting search to SMILES:', newState[0]);
              setSearch(newState[0]);
          }
          else if (representation === "SMARTS"){
              console.log('Setting search to SMARTS:', newState[1]);
              setSearch(newState[1]);
          }
          setToggleRepresentation(true);

          // Search will be triggered automatically by searchString dependency
        } else {
          console.log('Empty structure received from Ketcher');
        }
      };
      
    
    function switchRepresentations(event) {
        // Switch representations between SMARTS and smiles
      if (event === true) {
        setRepresentation("smiles");
        setSwitchCheck(true);
        setSearch(smiles);
      }
      // Remove label otherwise
      else {
        setRepresentation("SMARTS");
        setSwitchCheck(false);
        setSearch(SMARTS);
      }

      newSearch();
    }
    
    // loadmore
    function loadMore() {
        setSkip(skip => skip + interval);
        setSearchPage( searchPage => searchPage + 1);
        setIsLoadingMore(true);
        setSearchToggle(!searchToggle)
    }

    function newSearch() {
        setSkip(0);
        setSearchPage(1);
        setSVGResults([]);
        setResults([]);
        setIsLoading(true);
        
        // Force a search by triggering the search effect
        // Since searchString dependency handles this, we just need to ensure it fires
        const controller = new AbortController();
        const signal = controller.signal;
        
        setUpdatedParameters(false);
        loadImages(signal);
    }
 
    function loadImages(signal) {

        const fetchData = async () => {
            const molecule_data = await substructureSearch(searchString, interval, skip, signal);
            const svg_data = await retrieveAllSVGs(molecule_data, searchString, signal);

            return [ molecule_data, svg_data ]
        }

        fetchData()
        .then( (items )=> {
            if (searchPage == 1) {
            setSVGResults(items[1]);
            setResults(items[0]);
            }

            else {
                setSVGResults(svg_results.concat(items[1]));
                setResults(results.concat(items[0]) )
            }

            setIsLoading(false);
            setIsLoadingMore(false);
            setValidSmiles(true);

          })
        .catch( (error) => {
            console.log(error)
            setValidSmiles(false);
            setResults([]);
            setSVGResults([])
            setIsLoading(false)
            setIsLoadingMore(false) 
        } )

    }

    // If search string changes, then parameters have been updated.
    useEffect(() => {
      setUpdatedParameters(true);
    }, [searchString])

    // initial load of data
    // and load when search changes. 
    useEffect( ( ) => {
        // Only run search if we have a search string
        if (!searchString) return;
        
        const controller = new AbortController();
        const signal = controller.signal;

        setUpdatedParameters(false);
        setIsLoading(true);
        loadImages(signal);
        
        return () => {
          controller.abort();
        }
      },
        // eslint-disable-next-line react-hooks/exhaustive-deps
        [ searchString ] 
    );

      const _handleKeyDown = useCallback(
        (event) => {
          if (event.key === "Enter") {
            setFromKetcher(false);
            newSearch();
          }
        },
        [newSearch, setFromKetcher]
      );

      const searchButton = useCallback(() => {
        setFromKetcher(false);
        newSearch();
      }, [newSearch, setFromKetcher]);

    return (
        <Container maxWidth="lg">
        <Typography variant="h2" textAlign="center">Substructure Search</Typography>

        <Box display="flex" justifyContent="center" pt={2}>
        <FullScreenDialog ketcherCallBack={ketcherCallBack} />
        </Box>

        <Box display="flex" justifyContent="center">
        <TextField id="search-outline" 
                  label="Search for SMILES or SMARTS String" 
                  variant="outlined"
                  value= {searchString} 
                  onChange = { event => setSearch( event.target.value ) }
                  onKeyDown = { (e) => _handleKeyDown(e) }
                  InputProps={{endAdornment: <Button onClick={ () => { searchButton() } } 
                  >
                    Search
                    </Button>}}
                    />
            
        </Box>

        { toggleRepresentation &&
        <Grid component="label" container alignItems="center" spacing={1} sx={{position: 'flex', flexDirection: 'row', justifyContent: 'center', alignItems:'center'}}>
            <Grid item>SMARTS</Grid>
            <Grid item>
                <Switch
                checked={ switchCheck }
                onChange={ event => switchRepresentations(event.target.checked)}
                disabled={!fromKetcher}
                />
            </Grid>
            <Grid item>SMILES</Grid>
        </Grid>
        }
      <Box sx={{ display: 'flex', justifyContent: 'center', my: 3 }}>
                  { isLoadingMore ? <CircularProgress /> 
                  : 
                    <>
                      <Button disabled={updatedParameters} variant="contained" sx={{ my: 3 }} onClick={ () => loadMore() } >
                      <span style={{ textTransform: 'capitalize', fontSize: '16px' }}>
                        Load More
                      </span>
                      </Button>
                      <DropDownButton molecule_ids={moleculeIDs} dataTypes={availableDataTypes} />
                    </>
                  }          
      </Box>
        <Container sx={{display: 'flex', justifyContent: 'center', my: 3}}>
            <Box sx={{ display: 'flex' }}>
             { isLoading && <CircularProgress /> }
             { !isLoading && !validSmiles  && <Typography>Search string is not valid SMILES or SMARTS. Please provide valid SMILES or SMARTS strings.</Typography> }
             { !isLoading && validSmiles && Object.keys(svg_results).length > 0 && 
             <Container>
              <Container sx={{display: 'flex', justifyContent: 'left', my: 3}}>    
                  <Typography  sx={{ fontStyle: 'italic' }}>Showing {interval + interval * (searchPage-1)} results.</Typography>
              </Container> 
                { dynamicGrid(svg_results)  }
                <Box sx={{ display: 'flex', justifyContent: 'center', my: 3 }}>
                  { isLoadingMore ? <CircularProgress /> 
                  : 
                    <>
                      <Button disabled={updatedParameters} variant="contained" sx={{ my: 3 }} onClick={ () => loadMore() } >
                      <span style={{ textTransform: 'capitalize', fontSize: '16px' }}>
                        Load More
                      </span>
                      </Button>
                      <DropDownButton molecule_ids={moleculeIDs} dataTypes={availableDataTypes} />
                    </>
                  }          
                </Box>
               </Container>  }
            { !isLoading && validSmiles && Object.keys(svg_results).length==0 && <Typography>No results found for SMILES string.</Typography> } 
            </Box>
        </Container>

        </Container>
    )
}
