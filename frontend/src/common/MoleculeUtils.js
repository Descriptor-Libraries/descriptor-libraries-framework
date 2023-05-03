import React from 'react';
import { styled } from '@mui/material/styles';
import { Typography } from "@mui/material";
import Paper from '@mui/material/Paper';
import Grid from '@mui/material/Grid';
import Container from '@mui/material/Container';

import { createTheme } from '@mui/material/styles';

const theme = createTheme({
  palette: {
    primary: {
      main: "#ed1c24",
    }
  },
});

const Item = styled(Paper)(({ theme }) => ({
    backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
    ...theme.typography.body2,
    padding: theme.spacing(1),
    textAlign: 'center',
    color: theme.palette.text.secondary,
  }));

async function substructureSearch(substructure, limit=48, skip=0, signal) {
    let encoded = encodeURIComponent(substructure);
    
    const response =  await fetch(`/api/molecules/search/?substructure=${encoded}&skip=${skip}&limit=${limit}`, {signal: signal})

    if (!response.ok) {
        throw new Error('invalid smiles')
    }

    else {
        return await response.json()
    }
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
        <Grid container spacing={2} sx= {{ mt: 6 }}>
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
            
            {result.distance !== undefined && (
                <Typography sx={{ wordBreak: "break-word" }}>
                  <strong>Distance: </strong> {result.distance.toExponential(2)}
                </Typography>
              )}

            </Item>} 
        </Grid>
        ))
        }
        
    </Grid>
    </Container>
    )
  }


async function retrieveSVG(smiles, substructure = undefined, distance = undefined, signal) {
    /**
     * Retrieve SVG representation of a molecule.
     *
     * @param {string} smiles - SMILES representation of the molecule.
     * @param {string} [substructure=undefined] - SMILES representation of the substructure to highlight (optional).
     * @param {number} [distance=undefined] - Distance between the molecule and a reference molecule (optional).
     * @param {AbortSignal} [signal] - Signal object for aborting the fetch request (optional).
     * @returns {Promise<Object>} - A Promise that resolves to an object containing the SVG, SMILES, and distance (if provided).
     */
    let encoded = encodeURIComponent(smiles);
    let encodedSub = encodeURIComponent(substructure);
  
    const response = await fetch(`/depict/cow/svg?smi=${encoded}&sma=${encodedSub}`, {signal: signal});
  
    const svg = await response.text();
    let result = {}
    result["svg"] = svg;
    result["smiles"] = smiles;
  
    if (typeof distance !== "undefined") {
      result["distance"] = distance;
    }
  
    return result;
  }
  

  async function retrieveAllSVGs(items, substructure = undefined) {
    /**
   * Retrieve SVGs for a list of molecules.
   *
   * @param {Array<Object>} items - An array of objects containing molecule data, including SMILES and distance (if applicable).
   * @param {string} [substructure=undefined] - SMILES representation of the substructure to highlight in all molecules (optional).
   * @returns {Promise<Array<Object>>} - A Promise that resolves to an array of objects containing SVG, SMILES, and distance for each molecule.
   */
    return await Promise.all(items.map((item) => {
      return retrieveSVG(item.smiles, substructure, item.dist);
    }));
  }
  
  


export { retrieveAllSVGs, dynamicGrid, substructureSearch, theme };