import Graph from "../components/Graph"
import React, { useEffect, useState, useMemo } from 'react';
import { useParams } from "react-router-dom";
import { Box, Grid, Container, TextField, MenuItem, Card, CardContent, Select, InputLabel, FormControl} from "@mui/material";
import { DataGrid, GridFooterContainer, GridFooter } from "@mui/x-data-grid";
import Typography from '@mui/material/Typography';
import { CircularProgress } from "@mui/material";
import { retrieveSVG } from "../common/MoleculeUtils";
import { Stage, Component } from "react-ngl";

export default function MoleculeInfo() {
   const params = useParams();
   const [ molData, setMolData ] = useState([]);
   const [ umapNeighborData, setUmapNeighborData ] = useState([]);
   const [ pcaNeighborData, setPcaNeighborData ] = useState([]);
   const [ identifierData, setIdentifierData ] = useState([]);
   const [ neighborData, setNeighborData ] = useState([]);
   const [ components, setComponents ] = useState(["1", "2"]);
   const [ type, setType ] = useState("umap");
   const [ svg, setSvg ] = useState({});
   const [ allConformers, setAllConformers ] = useState([]);
   const [ conformer, setConformer ] = useState("");

   const reprList = useMemo(() => [{
      type: 'ball+stick'
    }], []);

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

   function switchDimensionality(event) {
      setType(event.target.value);
      console.log(event.target.value)
      if (event.target.value === "umap")
      {
         setComponents(["1", "2"]);
         setNeighborData(umapNeighborData);
      }
      else {
         setComponents(["1", "2", "3", "4"]);
         setNeighborData(pcaNeighborData);
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
            const umap_neighbor_data = await dimensionality(molid, "umap", ["1", "2"], signal);
            const pca_neighbor_data = await dimensionality(molid, "pca", ["1", "2", "3", "4"], signal);
            const svg_data = await retrieveSVG(molecule_data.smiles, signal);
            const identifier_data = await identifiers(molecule_data.smiles, signal);
            return [ molecule_data, umap_neighbor_data, pca_neighbor_data, svg_data, identifier_data ]
         }

         fetchData()
         .catch( (error) => {
            console.log(error);
         })
         .then( (items )=> {
            setMolData(items[0]);
            setUmapNeighborData(items[1]);
            setPcaNeighborData(items[2]);
            setSvg(items[3]);
            setIdentifierData(items[4][0]);
            // Initial set neighbor data to umap so it appears on load.
            setNeighborData(items[1]);
            // If we have conformers, we can set our states.
            if (items[0].conformers_id.length > 0)
            {
               // Get list of conformers
               setAllConformers(items[0].conformers_id);
               // Set conformer to the first one available.
               setConformer(items[0].conformers_id[0].toString());
            }
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

   useEffect( ( ) => {
      console.log(conformer);
   },
      [ conformer ]
   );

   return (
      <Container maxWidth="xl" sx={{ display: 'flex', alignItems: 'center' }}>
         <Grid container spacing={2} maxWidth="xl" sx={{alignItems: 'center'}}>
            <Grid item xs={6} sx={{mt: 3}}>
                  {Object.keys(svg).length > 0 && <Box component="img" alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(svg.svg)}`}></Box>}
                  {Object.keys(molData).length > 0 && 
                        <Card>
                           <CardContent>
                           <Typography align='left'> <strong>Smiles:</strong> {molData.smiles} </Typography>
                           <Typography align='left'> <strong>InChI:</strong> {identifierData.InChI} </Typography>
                           <Typography align='left'> <strong>InChIKey:</strong> {identifierData.InChIKey} </Typography>
                           <Typography align='left'> <strong>Molecular Weight:</strong> {molData.molecular_weight.toFixed(2)} </Typography>
                           </CardContent>
                        </Card>}
            </Grid>
            <Grid item xs={6}>
               {Object.keys(molData).length > 0 && Table(molData.ml_data)}
            </Grid>
            <Grid item xs={allConformers.length > 0 && conformer.length > 0 ? 6 : 0}>
               {allConformers.length > 0 && conformer.length > 0 && 
               <Container>
                  <FormControl fullWidth variant="standard">
                     <InputLabel id="conformer">Conformer</InputLabel>
                     <Select
                     labelId="conformer"
                     id="dimension-outline"
                     value={conformer}
                     style={{width: 125}}
                     MenuProps={{ PaperProps: { sx: { maxHeight: 200 } } }}
                     onChange={ function(event) {setConformer(event.target.value.toString());} }
                     >
                        {allConformers.map((item, index) => (
                           <MenuItem key={index} value={item}>{item}</MenuItem>
                        ))}
                     </Select>
                  </FormControl>
                  <Box
                  display="flex"
                  justifyContent="center"
                  alignItems="center"
                  >
                     <Stage width="600px" height="600px" params={{backgroundColor: 'white'}} cameraState={{distance: -20}}>
                        <Component path={"/api/conformers/export/"+ conformer + ".sdf"} reprList={reprList} />
                     </Stage>
                  </Box>
               </Container>}
            </Grid>
            <Grid item xs={allConformers.length > 0 && conformer.length > 0 ? 6 : 12}>
                  {Object.keys(neighborData).length > 0 ? (<TextField
                  sx={{ mb: 0.5 }}
                  select
                  id="dimension-outline"
                  value={type}
                  onChange={event => switchDimensionality(event)}
                  >
                        <MenuItem value={"umap"}>UMAP</MenuItem>
                        <MenuItem value={"pca"}>PCA</MenuItem>
                  </TextField>)
                  :
                  <CircularProgress sx={{ color: "#ed1c24" }} />
                  }
                  {Object.keys(neighborData).length > 0 && 
                  <Container sx={{ display: 'flex', height: 600, mb: 10}}>
                     <Graph molData={neighborData} componentArray={components} type={type} neighborSearch={true}></Graph>
                  </Container>}
            </Grid>
         </Grid>
      </Container>
   )
}
