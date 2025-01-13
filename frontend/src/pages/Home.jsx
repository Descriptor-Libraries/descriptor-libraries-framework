import React, { useEffect, useState } from 'react';

import { Typography, Box, Grid } from "@mui/material";
import { useTheme } from '@mui/material/styles';

import Graphv2 from '../components/Graphv2';  

import SearchIcon from '@mui/icons-material/Search';
import BubbleChartIcon from '@mui/icons-material/BubbleChart';
import DownloadIcon from '@mui/icons-material/Download';
//import AutoStoriesIcon from '@mui/icons-material/AutoStories';
import InfoIcon from '@mui/icons-material/Info';
import StatsGrid from '../components/StatsGrid';
import IconLink from '../components/IconLink';

import CheckFileAndDownloadIcon from '../components/FileDownload';

import purify from 'dompurify';

const Badge = ({ isMobile }) => {
  const displayStyle = {
    maxWidth: isMobile ? '120px' : '240px', // This sets the maximum width of the image
    height: '15vh' // Set the height of the image to 20% of the viewport height (box is 35% of vh)
  };

  return (
    <img src={`/${document.location.pathname.split('/')[1]}/brand/logo.svg`} style={displayStyle} alt="logo" />
  );
};


function Home() {
   const [ molData, setMolData ] = useState([]);
   const [isMobile, setIsMobile] = useState(window.innerWidth < 768);
   const [name, setName] = useState("");
   const [tagline, setTagline] = useState("");
   const [visType, setVisType] = useState(""); 
   const theme = useTheme();

   useEffect(() => {
    async function initialize() {

      // First get the configuration
      const response = await fetch(`/${document.location.pathname.split('/')[1]}/brand/names.json`);
      const data = await response.json();
      
      // Set the basic info
      setName(data[0].name);
      setTagline(data[0].tagline);
      setVisType(data[0]?.default_display ? 
        (data[0].default_display === 'pca' ? 'pc' : data[0].default_display) 
        : 'umap'
      );

      // Then get the molecule data
      const molecule_data = await dimensions(data[0]?.default_display ?? 'umap');
      setMolData(molecule_data);
    }
  
    initialize().catch(error => console.error(error));
  }, []); // Run once on mount

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


   async function dimensions(componentType) {
      /**
       * Requests general umap data from the backend.
       * @return {json}  The response json.
       */
         let encoded = encodeURIComponent("1,2");

         const response =  await fetch(`/api/${document.location.pathname.split('/')[1]}/molecules/dimensions/?type=${componentType}&components=${encoded}&limit=10000`)
      
         if (!response.ok) {
            throw new Error('Invalid Molecule Id')
         }
      
         else {
            return await response.json()
         }
   }

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
        <Badge isMobile={isMobile} />
        <Typography variant={ isMobile ? "h4" :"h2"} color="white">
            { name }
        </Typography>
    </Box>
      {
        <Typography variant={ isMobile ? "subtitle1" :"h5"} color="white" textAlign="center" dangerouslySetInnerHTML={{ __html:purify.sanitize(tagline) }} >
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
          <CheckFileAndDownloadIcon />
        {/* Hide until we're ready to add
          <Grid item>
            <IconLink IconElement={AutoStoriesIcon} text="Documentation" link="/docs" reloadDocument></IconLink>
          </Grid>
          */}
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
  {!isMobile && <Graphv2 molData={molData} componentType={visType} containerStyle={{ width: '100%', height: '90%' }}></Graphv2>}
</Box>
    </>
  );
}

export default Home;
