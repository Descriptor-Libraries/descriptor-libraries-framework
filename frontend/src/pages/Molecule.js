import Graph from "../components/Graph"
import React, { useEffect, useState } from 'react';
import { useParams } from "react-router-dom";
import { Box, Grid, Item } from "@mui/material";
import Typography from '@mui/material/Typography';
import { retrieveSVG } from "../common/MoleculeUtils";

export default function MoleculeInfo() {
   const params = useParams();
   const [ molData, setMolData ] = useState([]);
   const [ neighborData, setNeighborData ] = useState([]);
   const [ components, setComponents ] = useState(["1", "2"]);
   const [ type, setType ] = useState("umap");
   const [ svg, setSvg ] = useState({});

   async function molecule(molecule_id, signal) {
      /**
       * Requests general umap or pca data from the backend.
       * @param {number} molecule_id Id of the molecule to search on.
       * @param {AbortSignal} signal Abortsignal object.
       */
         const response =  await fetch(`/api/molecules/${molecule_id}`, {signal: signal})
      
         if (!response.ok) {
            throw new Error('Invalid Molecule Id')
         }
      
         else {
            return await response.json()
         }
   }

   async function dimensionality(molecule_id, type, components, signal, limit=10) {
      /**
       * Requests general umap or pca data from the backend.
       * @param {number} molecule_id Id of the molecule to search on.
       * @param {string} type Type of dimensionality reduction. Can be one of PCA or UMAP.
       * @param {string} components String of comma separated integers.
       * @param {AbortSignal} signal Abortsignal object.
       * @param {number} limit Limit of the search.
       * @return {json}  The response json.
       */
         let encoded = encodeURIComponent(components);

         const response =  await fetch(`/api/molecules/${molecule_id}/neighbors/?type=${type}&components=${encoded}&skip=0&limit=${limit}`, {signal: signal})
      
         if (!response.ok) {
            throw new Error('Invalid Molecule Id')
         }
      
         else {
            return await response.json()
         }
   }

   function loadData(signal, molid) {
      /**
       * Main driver function which loads the neighbors for a molecule requested by the user.
       * @param {AbortSignal} signal Abortsignal object.
       */
         const fetchData = async () => {
            const molecule_data = await molecule(molid, signal);
            const neighbor_data = await dimensionality(molid, type, components, signal);
            const svg_data = await retrieveSVG(molecule_data.smiles, signal);
            return [ molecule_data, neighbor_data, svg_data ]
         }

         fetchData()
         .catch( (error) => {
            console.log(error);
         })
         .then( (items )=> {
            setMolData(items[0]);
            setNeighborData(items[1]);
            setSvg(items[2]);
      })
      }
   
      // initial load of data
      // and load when search changes. 
      useEffect( ( ) => {
         const controller = new AbortController();
         const signal = controller.signal;

         // setUpdatedParameters(false);
         loadData(signal, params.molid);

         return () => {
         controller.abort();
         }
      },
         // eslint-disable-next-line react-hooks/exhaustive-deps
         [ params ]
      );

   return (
      <Grid container rowSpacing={1}>
         <Grid item xs={6}>
               {Object.keys(svg).length > 0 && <Box sx={{ my: 3 }} component="img" alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(svg.svg)}`}></Box>}
               {Object.keys(molData).length > 0 && <Box sx={{ my: 3 }}>
                  <Typography> Smiles: {molData.smiles} </Typography>
                  <Typography> Molecular Weight: {molData.molecular_weight.toFixed(2)} </Typography>
               </Box>}
         </Grid>
         <Grid item xs={6}>
            Table of Information
         </Grid>
         <Grid item xs={6}>
            Conformer Structures
         </Grid>
         <Grid item xs={6}>
               {Object.keys(neighborData).length > 0 && <Graph molData={neighborData} componentArray={components} type={type} neighborSearch={true}></Graph>}
         </Grid>
      </Grid>
      // <Container maxWidth="xl" sx={{display: 'flex', flexDirection: "column", height: 850, alignItems: 'center'}}>
         
         
      // </Container>
   )
}
