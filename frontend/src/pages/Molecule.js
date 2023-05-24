import Graph from "../components/Graph"
import React, { useEffect, useState } from 'react';
import { Container } from "@mui/material";
import { Typography } from "@mui/material";

function Molecule() {
   const [ molData, setMolData ] = useState([]);
   const [ components, setComponents ] = useState(["1", "2"]);
   const [ type, setType ] = useState("umap");

   async function dimensionality(molecule_id, type, components, signal, limit=10) {
      /**
       * Requests general umap or pca data from the backend.
       * @param {string} type Type of dimensionality reduction. Can be one of PCA or UMAP.
       * @param {string} components String of comma separated integers.
       * @param {AbortSignal} signal Abortsignal object.
       * @param {number} limit Limit of the search.
       * @return {json}  The response json.
       */
         let encoded = encodeURIComponent(components);

         const response =  await fetch(`/api/molecules/${molecule_id}/neighbors/?type=${type}&components=${encoded}&limit=${limit}`, {signal: signal})
      
         if (!response.ok) {
            throw new Error('Invalid Molecule Id')
         }
      
         else {
            return await response.json()
         }
   }

   function loadData() {
      /**
       * Main driver function which loads the data.
       * 
       */
         const fetchData = async () => {
            const molecule_data = await dimensionality(type, components);
            return molecule_data
         }

         fetchData()
         .catch( (error) => {
            console.log(error) 
         })
         .then( (items )=> {   
            setMolData(items);   
      })
      }
   
      useEffect( ( ) => { 
      loadData() },  
      [ type ] 
   );

   return (
      <Container maxWidth="xl" sx={{display: 'flex', flexDirection: "column", height: 850, alignItems: 'center'}}>
         <Typography variant="h2">
         Kraken Webapp
         </Typography>
         {<Graph molData={molData} componentArray={components} type={type} neighborSearch={false}></Graph>}
      </Container>
   )
}

export default Molecule;
