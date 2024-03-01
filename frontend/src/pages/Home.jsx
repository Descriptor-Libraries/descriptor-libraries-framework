import React, { useEffect, useState } from 'react';

import { Typography, Box, Grid } from "@mui/material";
import { useTheme } from '@mui/material/styles';

import Graph from "../components/Graph";

import SearchIcon from '@mui/icons-material/Search';
import BubbleChartIcon from '@mui/icons-material/BubbleChart';
//import DownloadIcon from '@mui/icons-material/Download';
import AutoStoriesIcon from '@mui/icons-material/AutoStories';
import InfoIcon from '@mui/icons-material/Info';
import StatsGrid from '../components/StatsGrid';
import IconLink from '../components/IconLink';

import purify from 'dompurify';

const Badge = ({ isMobile }) => {
  const displayStyle = {
    maxWidth: isMobile ? '120px' : '240px', // This sets the maximum width of the image
    height: '15vh' // Set the height of the image to 20% of the viewport height (box is 35% of vh)
  };

  return (
    <img src={`${import.meta.env.BASE_URL}/brand/logo.svg`} style={displayStyle} alt="logo" />
  );
};


function Home() {
   const [ molData, setMolData ] = useState([]);
   const [isMobile, setIsMobile] = useState(window.innerWidth < 768);
   const [name, setName] = useState("");
   const [tagline, setTagline] = useState("");
   const theme = useTheme();

   useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}/brand/names.json`)
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.json();
        })
        .then(data => {
          setName(data[0].name);
          setTagline(data[0].tagline);
        })
        .catch(error => {
            console.error('Error fetching stats:', error);
        });
}, []);

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


   async function umap() {
      /**
       * Requests general umap data from the backend.
       * @return {json}  The response json.
       */
         let encoded = encodeURIComponent("1,2");

         const response =  await fetch(`/api${import.meta.env.BASE_URL}/molecules/dimensions/?type=umap&components=${encoded}&limit=1000`)
      
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
            const molecule_data = await umap();
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
   }, []);

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
  {!isMobile && <Graph molData={molData} componentArray={["1", "2"]} neighborSearch={false} containerStyle={{ width: '100%', height: '90%' }}></Graph>}
</Box>
    </>
  );
}

export default Home;
