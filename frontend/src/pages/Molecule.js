import Graph from "../components/Graph"
import React, { useEffect, useState } from 'react';
import { useParams } from "react-router-dom";
import { Box, Grid, Container } from "@mui/material";
import { DataGrid, GridFooterContainer, GridFooter } from "@mui/x-data-grid";
import Typography from '@mui/material/Typography';
import { retrieveSVG } from "../common/MoleculeUtils";

export default function MoleculeInfo() {
   const params = useParams();
   const [ molData, setMolData ] = useState([]);
   const [ neighborData, setNeighborData ] = useState([]);
   const [ identifierData, setIdentifierData ] = useState([]);
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

   async function identifiers(smiles, signal) {
      /**
       * Requests general umap or pca data from the backend.
       * @param {string} smiles Smiles of the molecule.
       * @param {AbortSignal} signal Abortsignal object.
       * @return {json}  The response json.
       */
         let encoded = encodeURIComponent(smiles);

         const response =  await fetch(`/api/molecules/identifiers/?smiles=${encoded}`, {signal: signal})
      
         if (!response.ok) {
            throw new Error('Invalid Molecule Smiles')
         }
      
         else {
            return await response.json()
         }
   }

   function Table(data) {
      let columns = [];
      let rows = [];

      let keys = Object.keys(data);
      let properties = Object.keys(data["max_data"])

      columns.push({field: "id", flex: 1, headerClassName: "super-app-theme--header"});
      
      // Loop through all the keys and create columns and rows. Avoid boltzmann_averaged_data since it does not have the same keys as the rest.
      for (const element of keys) {
         if (element != "boltzmann_averaged_data")
         {
            columns.push({field: element, flex: 0.75, headerClassName: "super-app-theme--header"});
         }
      }

      // Make the rows of the table
      for (const property of properties) {
         let newObj = {id: property};
         for (const element of keys) {
            if (element != "boltzmann_averaged_data")
            {
               newObj[element] = data[element][property];
            }
         }
         rows.push(newObj);
      }

      function CustomFooter () {
         return (
           <GridFooterContainer>
             <Typography sx= {{mx: 1, color: 'gray'}}> ML Data </Typography>
             <GridFooter sx={{
               border: 'none', // To delete double border.
               }} />
           </GridFooterContainer>
         );
       }

      return (
         <Box
            sx={{
            width: '100%',
            '& .super-app-theme--header': {
               backgroundColor: 'rgba(237, 28, 36)',
               color: 'white',
            },
            }}
         >
            <DataGrid
               disableColumnMenu
               rows={rows}
               columns={columns}
               components={{Footer: CustomFooter}}
               initialState={{
                  pagination: {
                     paginationModel: { page: 0, pageSize: 4 },
                  },
               }}
            />
         </Box>
      )
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
            const identifier_data = await identifiers(molecule_data.smiles, signal);
            return [ molecule_data, neighbor_data, svg_data, identifier_data ]
         }

         fetchData()
         .catch( (error) => {
            console.log(error);
         })
         .then( (items )=> {
            setMolData(items[0]);
            setNeighborData(items[1]);
            setSvg(items[2]);
            setIdentifierData(items[3][0]);
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
      <Container maxWidth="xl" sx={{ my: 5, display: 'flex', height: 850, alignItems: 'center'}}>
         <Grid container rowSpacing={1} maxWidth="xl" sx={{alignItems: 'center'}}>
            <Grid item xs={6}>
                  {Object.keys(svg).length > 0 && <Box sx={{ my: 3 }} component="img" alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(svg.svg)}`}></Box>}
                  {Object.keys(molData).length > 0 && <Box sx={{ my: 3 }}>
                     <Typography> <strong>Smiles:</strong> {molData.smiles} </Typography>
                     <Typography> <strong>InChI:</strong> {identifierData.InChI} </Typography>
                     <Typography> <strong>InChIKey:</strong> {identifierData.InChIKey} </Typography>
                     <Typography> <strong>Molecular Weight:</strong> {molData.molecular_weight.toFixed(2)} </Typography>
                  </Box>}
            </Grid>
            <Grid item xs={6}>
               {Object.keys(molData).length > 0 && Table(molData.ml_data)}
            </Grid>
            <Grid item xs={12}>
                  {Object.keys(neighborData).length > 0 && <Graph molData={neighborData} componentArray={components} type={type} neighborSearch={true}></Graph>}
            </Grid>
         </Grid>
      </Container>
   )
}
