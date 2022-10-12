import React, { useEffect, useState, useCallback } from 'react';
import { styled } from '@mui/material/styles';
import { TextField, Typography } from "@mui/material";
import Paper from '@mui/material/Paper';
import Grid from '@mui/material/Grid';
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import CircularProgress from '@mui/material/CircularProgress';

import Button from '@mui/material/Button';

const Item = styled(Paper)(({ theme }) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.text.secondary,
  }));

async function substructureSearch(substructure, limit=48, skip=0) {
    let encoded = encodeURIComponent(substructure);

    const response =  await fetch(`/api/v1/molecule/search/?substructure=${encoded}&skip=${skip}&limit=${limit}`)

    console.log(`/api/v1/molecule/search/?substructure=${encoded}&skip=${skip}&limit=${limit}`)

    if (!response.ok) {
        throw 'invalid smiles'
    }

    else {
        return await response.json()
    }
}

async function retrieveSVG( smiles, substructure ) {
    let encoded = encodeURIComponent(smiles);
    let encodedSub = encodeURIComponent(substructure);

    const response = await fetch(`/depict/cow/svg?smi=${encoded}&sma=${encodedSub}&zoom=1.25&w=50&h=50`);

    const svg = await response.text();
    let result = {}
    result["svg"] = svg;
    result["smiles"] = smiles;
    return result
}

async function retrieveAllSVGs( items, substructure ) {
    return await Promise.all( items.map( (item) => { 
        return retrieveSVG(item.smiles, substructure)
     } ) )
}

function dynamicGrid( svgs ) {

    return (
        <Container>
        <Grid container spacing={2} sx= {{ mt: 3 }}>
        {
        svgs.map((result) => (
        <Grid item xs={12} md={4}>
            <Item>
            <img alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(result.svg)}`} />
            <Typography>{ result.smiles }</Typography>
            </Item> 
        </Grid>
        ))
        }
        
    </Grid>
        <Button variant="contained" style={{backgroundColor: "#ed1c24"}} sx={{ my: 3 }}>Load More</Button>
    </Container>
    )
}

export default function SearchHook () {

    const interval = 21

    const [ searchString, setSearch ] = useState('PC=C');
    const [ skip, setSkip ] = useState(0);
    const [ limit, setLimit ] = useState(interval);
    const [ results, setResults ] = useState([]);
    const [ validSmiles, setValidSmiles ] = useState(true);
    const [ svg_results, setSVGResults ] = useState([])
    const [ searchPage, setSearchPage ] = useState(0);
    const [ isLoading, setIsLoading ] = useState(true);

    // 
    function loadImages() {
        setIsLoading(true);

        const fetchData = async () => {
            const molecule_data = await substructureSearch(searchString, limit, skip);
            const svg_data = await retrieveAllSVGs(molecule_data, searchString);

            return [ molecule_data, svg_data ]
        }

        fetchData()
        .catch( (error) => { 
            console.log(error)
            setValidSmiles(false);
            setResults([]); 
            setIsLoading(false);
        } )
        .then( (items )=> {
            setSVGResults(items[1]);
            setResults(items[0]);
            setIsLoading(false);
          })

    }

    // initial load of data and what happens when search
    // string changes
    useEffect( () => { loadImages() }, [] );

    function _handleKeyDown(event) {
        if (event.key === "Enter") {
          loadImages();
        }
      }

    return (
        <Container maxWidth="lg">
        <h2>Substructure Search</h2>
        <TextField id="search-outline" 
                  label="Enter a SMILES String to Search" 
                  variant="outlined"
                  defaultValue= {searchString} 
                  onChange = { event => setSearch( event.target.value ) }
                  onKeyDown = { (e) => _handleKeyDown(e) }
                  InputProps={{endAdornment: <Button onClick={ () => { loadImages() } } 
                  >
                    Search
                    </Button>}}
                    />

        <Container sx={{display: 'flex', justifyContent: 'center', my: 3}}>
            <Box sx={{ display: 'flex' }}>
             { isLoading && <CircularProgress sx={{ color: "#ed1c24" }} /> }
             { !isLoading && !validSmiles  && <Typography>Invalid Smiles String</Typography> }
             { !isLoading && validSmiles && Object.keys(svg_results).length && <Container> { dynamicGrid(svg_results) } </Container> } 
            </Box>
        </Container>

            

        </Container>
    )
}
