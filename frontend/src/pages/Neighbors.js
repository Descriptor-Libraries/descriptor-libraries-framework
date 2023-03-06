import React, { useEffect, useState } from 'react';
import { styled } from '@mui/material/styles';
import { TextField, Typography } from "@mui/material";
import Paper from '@mui/material/Paper';
import Grid from '@mui/material/Grid';
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import CircularProgress from '@mui/material/CircularProgress';

import Button from '@mui/material/Button';
import MenuItem from '@mui/material/MenuItem';

const Item = styled(Paper)(({ theme }) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.text.secondary,
  }));

async function NeighborSearch(molecule_id, type, components, limit=48, skip=0) {
    let encoded = encodeURIComponent(components);

    const response =  await fetch(`/api/molecules/${molecule_id}/neighbors/?type=${type}&components=${encoded}&skip=${skip}&limit=${limit}`)

    if (!response.ok) {
        throw new Error('Invalid Molecule Id')
    }

    else {
        return await response.json()
    }
}

export default function NeighborSearchHook () {

    const interval = 15;

    const [ moleculeid, setSearch ] = useState(1);
    const [ type, setType ] = useState("pca");
    const [ components, setComponents ] = useState("1,2,3,4");
    const [ skip, setSkip ] = useState(1);
    const [ results, setResults ] = useState([]);
    const [ validMolecule, setValidMolecule ] = useState(true);
    const [ svg_results, setSVGResults ] = useState([])
    const [ searchPage, setSearchPage ] = useState(1);
    const [ isLoading, setIsLoading ] = useState(true);
    const [ searchToggle, setSearchToggle ] = useState(true);
    const [ isLoadingMore, setIsLoadingMore ] = useState(false);
    
    // loadmore
    function loadMore() {
        setSkip(skip => skip + interval);
        setSearchPage( searchPage => searchPage + 1);
        setIsLoadingMore(true);
    }

    function newSearch() {
        setSkip(1);
        setSearchPage(1);
        setSVGResults([]);
        setResults([]);
        // Just need to toggle this to make sure it toggles
        // so that effect will be triggered
        setIsLoading(true);
        setSearchToggle(!searchToggle);
    }
 
    function loadImages() {

        const fetchData = async () => {
            const molecule_data = await NeighborSearch(moleculeid, type, components, interval, skip);
            return [ molecule_data ]
        }

        fetchData()
        .catch( (error) => {
            console.log(error) 
            setValidMolecule(false);
            setResults([]);
            setSVGResults([])
            setIsLoading(false)
            setIsLoadingMore(false) 
        } )
        .then( (items )=> {
            console.log(items);
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

          })

    }

    // initial load of data 
    useEffect( ( ) => { 
        loadImages() }, 

        // eslint-disable-next-line react-hooks/exhaustive-deps
        [ searchPage, searchToggle ] 
    );

    function _handleKeyDown(event) {
        if (event.key === "Enter") {
          loadImages();
        }
      }

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
                  onKeyDown = { (e) => _handleKeyDown(e) }
                  InputProps={{endAdornment: <Button onClick={ () => { newSearch() } } 
                  >
                    Search
                    </Button>}}
        />
        <TextField
            sx={{ m: 0.5 }}
            select
            id="dimension-outline"
            value={type}
            onChange={event => setType(event.target.value)}
        >
            <MenuItem value={"pca"}>PCA</MenuItem>
            <MenuItem value={"umap"}>UMAP</MenuItem>
        </TextField>
        <TextField
                  sx={{ m: 0.5}} 
                  id="component-outline" 
                  label="Enter Components" 
                  variant="outlined"
                  defaultValue= {components} 
                  onChange = { event => setComponents( event.target.value ) }
        />

        <Container sx={{display: 'flex', justifyContent: 'center', my: 3}}>
            <Box sx={{ display: 'flex' }}>
            { !isLoading && validMolecule && <Typography>No results found for SMILES string.</Typography> } 
            </Box>
        </Container>

            

        </Container>
    )
}
