import React, { useEffect, useState } from 'react';
import { styled } from '@mui/material/styles';
import { TextField, Typography } from "@mui/material";
import Paper from '@mui/material/Paper';
import Grid from '@mui/material/Grid';
import Box from '@mui/material/Box';
import Container from '@mui/material/Container';
import CircularProgress from '@mui/material/CircularProgress';

import Dialog from '@mui/material/Dialog';
import Slide from '@mui/material/Slide';

import { Ketcher } from 'ketcher-core';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Editor } from 'ketcher-react';
import "ketcher-react/dist/index.css";

import Button from '@mui/material/Button';

const structServiceProvider = new StandaloneStructServiceProvider();

const Transition = React.forwardRef(function Transition(props, ref) {
    return <Slide direction="up" ref={ref} {...props} />;
  });

const Item = styled(Paper)(({ theme }) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.text.secondary,
  }));

async function substructureSearch(substructure, limit=48, skip=0) {
    let encoded = encodeURIComponent(substructure);

    const response =  await fetch(`/api/molecules/search/?substructure=${encoded}&skip=${skip}&limit=${limit}`)

    if (!response.ok) {
        throw new Error('invalid smiles')
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
            <Typography sx={{ wordBreak: "break-word" }}>{ result.smiles }</Typography>
            </Item> 
        </Grid>
        ))
        }
        
    </Grid>
    </Container>
    )
}

export default function SearchHook () {

    const interval = 15;

    const [ searchString, setSearch ] = useState('PC=C');
    const [ skip, setSkip ] = useState(0);
    const [ results, setResults ] = useState([]);
    const [ validSmiles, setValidSmiles ] = useState(true);
    const [ svg_results, setSVGResults ] = useState([])
    const [ searchPage, setSearchPage ] = useState(1);
    const [ isLoading, setIsLoading ] = useState(true);
    const [ searchToggle, setSearchToggle ] = useState(true);
    const [ isLoadingMore, setIsLoadingMore ] = useState(false);
    const [ ketcher, setKetcher ] = useState();
    const [ smiles, setSmiles ] = useState();
    const [ SMARTS, setSMARTS ] = useState();

    // Ketcher
    function FullScreenDialog() {
        const [open, setOpen] = React.useState(false);
      
        const handleClickOpen = () => {
            setOpen(true);
        };
      
        const handleClose = () => {
            ketcher.getSmiles().then(result => {setSmiles(result);});
            setOpen(false);
        };
      
        return (
          <div>
            <Button onClick={handleClickOpen}>
              DRAW
            </Button>
            <Dialog
              fullWidth={true}
              maxWidth={"lg"}
              open={open}
              onClose={handleClose}
              TransitionComponent={Transition}
            >
                <Editor
                    staticResourcesUrl={process.env.PUBLIC_URL}
                    structServiceProvider={structServiceProvider}
                    onInit={(ketcher) => {
                        setKetcher(ketcher)
                    }}
                />
            </Dialog>
          </div>
        );
      }
    
    // loadmore
    function loadMore() {
        setSkip(skip => skip + interval);
        setSearchPage( searchPage => searchPage + 1);
        setIsLoadingMore(true);
    }

    function newSearch() {
        setSkip(0);
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
            const molecule_data = await substructureSearch(searchString, interval, skip);
            const svg_data = await retrieveAllSVGs(molecule_data, searchString);

            return [ molecule_data, svg_data ]
        }

        fetchData()
        .catch( (error) => {
            console.log(error) 
            setValidSmiles(false);
            setResults([]);
            setSVGResults([])
            setIsLoading(false)
            setIsLoadingMore(false) 
        } )
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

    }

    // initial load of data 
    useEffect( ( ) => { 
        loadImages() }, 

        // eslint-disable-next-line react-hooks/exhaustive-deps
        [ searchPage, searchToggle ] 
    );

    // Update searchString if smiles changes
    useEffect(() => {
        setSearch(smiles);
      }, [smiles]);

    // Update searchString if SMARTS changes
    useEffect(() => {
        setSearch(SMARTS);
      }, [SMARTS]);

    // New search if searchString changes
    useEffect(() => {
        loadImages();
      }, [searchString]);

    function _handleKeyDown(event) {
        if (event.key === "Enter") {
          loadImages();
        }
      }

    return (
        <Container maxWidth="lg">
        <h2>Substructure Search</h2>
        <TextField id="search-outline" 
                  label="Enter a SMILES or SMARTS String to Search" 
                  variant="outlined"
                  defaultValue = {searchString} 
                  value = {searchString}
                  onChange = { event => setSearch( event.target.value ) }
                  onKeyDown = { (e) => _handleKeyDown(e) }
                  InputProps={{endAdornment: FullScreenDialog()}}
                    />

        <Container sx={{display: 'flex', justifyContent: 'center', my: 3}}>
            <Box sx={{ display: 'flex' }}>
             { isLoading && <CircularProgress sx={{ color: "#ed1c24" }} /> }
             { !isLoading && !validSmiles  && <Typography>Search string is not valid SMILES or SMARTS. Please provide valid SMILES or SMARTS strings.</Typography> }
             { !isLoading && validSmiles && Object.keys(svg_results).length > 0 && 
             <Container> 
                { dynamicGrid(svg_results)  }
                { isLoadingMore ? <CircularProgress sx={{ color: "#ed1c24" }} /> : <Button variant="contained" style={{backgroundColor: "#ed1c24"}} sx={{ my: 3 }} onClick={ () => loadMore() }>Load More</Button> }
            </Container>  }
            { !isLoading && validSmiles && Object.keys(svg_results).length==0 && <Typography>No results found for SMILES string.</Typography> } 
            </Box>
        </Container>

            

        </Container>
    )
}
