import React, { useEffect, useState } from 'react';

import { Typography, Box, Grid } from "@mui/material";
import { useTheme } from '@mui/material/styles';

import OriginalKraken from '../common/OriginalKraken';
import Graph from "../components/Graph";

import SearchIcon from '@mui/icons-material/Search';
import BubbleChartIcon from '@mui/icons-material/BubbleChart';
//import DownloadIcon from '@mui/icons-material/Download';
import AutoStoriesIcon from '@mui/icons-material/AutoStories';
import InfoIcon from '@mui/icons-material/Info';
import StatsGrid from '../components/StatsGrid';
import IconLink from '../components/IconLink';


function Home() {
   const [ molData, setMolData ] = useState([]);
   const [ components, setComponents ] = useState(["1", "2"]);
   const [ type, setType ] = useState("umap");
   const [isMobile, setIsMobile] = useState(window.innerWidth < 768);
   const theme = useTheme();

   useEffect(() => {
    function checkMobile() {
      setIsMobile(window.innerWidth < 768);
    }

    // Set isMobile at the start in case it's not the initial render
    checkMobile();

    window.addEventListener('resize', checkMobile);

    // Cleanup the listener when the component is unmounted
    return () => window.removeEventListener('resize', checkMobile);
  }, []); // Empty array means this effect runs once on mount and cleanup on unmount


   async function umap(type, components, limit=1000) {
      /**
       * Requests general umap data from the backend.
       * @param {string} type Type of dimensionality reduction. Can be one of PCA or UMAP.
       * @param {string} components String of comma separated integers.
       * @param {number} limit Limit of the search.
       * @return {json}  The response json.
       */
         let encoded = encodeURIComponent(components);

         const response =  await fetch(`/api/molecules/dimensions/?type=${type}&components=${encoded}&limit=${limit}`)
      
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
            const molecule_data = await umap(type, components);
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
   useEffect(() => {
      loadData()
   }, [type]);

  return (
  <>
   <Box sx={{ backgroundColor: theme.palette.primary.main, 
            display: 'flex', 
            flexDirection: 'column', 
            alignItems: 'center', 
            justifyContent: 'space-around', 
            padding: 2, height: !isMobile && '35vh' }}>  
   <Box 
        sx={{ 
            display: 'flex', 
            alignItems: 'center', 
            gap: 2,
            flexDirection: isMobile ? 'column' : 'row',
        }}
    >
        <OriginalKraken sx={{ color: 'white', fontSize: isMobile ? '120px' : '160px' }} />
        <Typography variant={ isMobile ? "h4" :"h2"} color="white">
            KRAKEN
        </Typography>
    </Box>
      {
        <Typography variant={ isMobile ? "subtitle1" :"h5"} color="white" textAlign="center">
          <b>K</b>olossal vi<b>R</b>tual d<b>A</b>tabase
                  for mole<b>K</b>ular d<b>E</b>scriptors of orga<b>N</b>ophosphorus
                  ligands.
        </Typography>
      }
       <Grid container alignItems="center" justifyContent="center" spacing={3} sx={{ mb: 1 }}>
          <Grid item>
            <IconLink IconElement={SearchIcon} text="Substructure Search" link='/search'></IconLink>
          </Grid>
          <Grid item>
            <IconLink IconElement={BubbleChartIcon} text="Neighbor Search" link='/neighbors'></IconLink>
          </Grid>
          <Grid item>
            <IconLink IconElement={InfoIcon} text="Library Details" link="/library_details"></IconLink>
          </Grid>
          {/* Hide until we're ready to add
            <Grid item>
              <IconLink IconElement={DownloadIcon} text="Download" link="/download"></IconLink>
            </Grid> 
            */}
          <Grid item>
            <IconLink IconElement={AutoStoriesIcon} text="Documentation" link="/docs" reloadDocument></IconLink>
          </Grid>
        </Grid>
     </Box>
     <StatsGrid />
      <Box sx={{ 
        width: '100%',
        height: '80vh', 
        mt: 3,
        alignItems: 'center',  
        justifyContent: 'center', 
      }}>
  {!isMobile && <Graph molData={molData} componentArray={components} type={type} neighborSearch={false} containerStyle={{ width: '100%', height: '90%' }}></Graph>}
</Box>
    </>
  );
}

export default Home;
