import React, { useEffect, useState } from 'react';
import { Link } from "react-router-dom";

import { Container, Typography, Box, Grid, IconButton } from "@mui/material";

import OriginalKraken from '../common/OriginalKraken';
import Graph from "../components/Graph";

import SearchIcon from '@mui/icons-material/Search';
import BubbleChartIcon from '@mui/icons-material/BubbleChart';
import DownloadIcon from '@mui/icons-material/Download';
import AutoStoriesIcon from '@mui/icons-material/AutoStories';
import StatCard from '../components/SummaryCard';

const stats = [
  {
    "number": "1,558",
    "description": "DFT Calculated Ligands"
  },

  {
    "number": "300,000+",
    "description": "ML calculated Ligands"
  },

  {
    "number": "576",
    "description": "unique substituents"
  },
  
  {
    "number": "190",
    "description": "DFT-level descriptors"
  },

  {
    "number": "21,437",
    "description": "unique molecules"
  },

  {
    "number": "13.8",
    "description": "average conformers per ligand"
  }, 
  
]

function Home() {
   const [ molData, setMolData ] = useState([]);
   const [ components, setComponents ] = useState(["1", "2"]);
   const [ type, setType ] = useState("umap");
   const [isMobile, setIsMobile] = useState(window.innerWidth < 768);

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

      useEffect( ( ) => { 
      loadData() },  
      [ type ] 
   );
   useEffect(() => {
      loadData()
   }, [type]);

  return (
   <>
   <Box sx={{ backgroundColor: "#393536", display: 'flex', flexDirection: 'column', alignItems: 'center', justifyContent: 'space-around', padding: 2, height: !isMobile && '35vh' }}>  
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
       <Grid container justifyContent="center" spacing={3} sx={{ mb: 1 }}>
          <Grid item>
            <Link to="/search">
              <IconButton color="inherit">
                <SearchIcon sx={{ color: 'white', fontSize: 40 }} />
              </IconButton>
            </Link>
            <Typography color="white" textAlign="center">
              Substructure Search
            </Typography>
          </Grid>
          <Grid item>
            <Link to="/neighbors">
              <IconButton color="inherit">
                <BubbleChartIcon sx={{ color: 'white', fontSize: 40 }} />
              </IconButton>
            </Link>
            <Typography color="white" textAlign="center">
              Neighbor Search
            </Typography>
          </Grid>
          <Grid item>
            <Link to="/download">
              <IconButton color="inherit">
                <DownloadIcon sx={{ color: 'white', fontSize: 40 }} />
              </IconButton>
            </Link>
            <Typography color="white" textAlign="center">
              Download
            </Typography>
          </Grid>
          <Grid item>
            <Link to="/docs/" reloadDocument>
              <IconButton color="inherit">
                <AutoStoriesIcon sx={{ color: 'white', fontSize: 40 }} />
              </IconButton>
            </Link>
            <Typography color="white" textAlign="center">
              Documentation
            </Typography>
          </Grid>
        </Grid>
     </Box>

    <Container maxWidth="xl" sx={{ display: "flex", flexDirection: "column", alignItems: "center" }}>
      <Grid container spacing={2} sx={{ mt: 3 }}>
            {stats.map((stat, index) => (
              <StatCard key={index} number={stat.number} caption={stat.description} size={150} />
            ))}
      </Grid>
      <Box sx={{ width: '100%', mt: 3 }}>
        {!isMobile && <Graph molData={molData} componentArray={components} type={type} neighborSearch={false}></Graph>}
      </Box>
    </Container>

    </>
  );
}

export default Home;
