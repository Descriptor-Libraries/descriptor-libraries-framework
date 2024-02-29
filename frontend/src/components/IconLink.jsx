import React from 'react';
import { Link } from 'react-router-dom';
import { Typography, Box, IconButton } from "@mui/material";
import SearchIcon from '@mui/icons-material/Search';

function IconLink({ IconElement, text, link, reloadDocument = false }) {
  const linkProps = reloadDocument ? { to: link, reloadDocument } : { to: link };

  return (
    <Link style={{ textDecoration: "none" }} {...linkProps}>
      <Box textAlign="center">
        <IconButton color="inherit" size="large">
          {IconElement ? <IconElement sx={{ color: 'white', fontSize: "5vh" }}/> : <SearchIcon sx={{ color: 'white', fontSize: "5vh" }}/>}
        </IconButton>
        <Typography color="white" textAlign="center">
          {text}
        </Typography>
      </Box>
    </Link>
  )
}

export default IconLink;
