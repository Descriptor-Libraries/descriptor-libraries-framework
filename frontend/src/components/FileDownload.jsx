import React, { useState, useEffect } from 'react';
import Grid from '@mui/material/Grid';
import DownloadIcon from '@mui/icons-material/Download';
import IconLink from './IconLink'; 

const CheckFileAndDownloadIcon = () => {
  const [fileExists, setFileExists] = useState(false);

  useEffect(() => {
    const checkFileExists = async () => {
      try {
        // Attempt to fetch the file
        const response = await fetch(`${document.location.pathname.split('/')[1]}/data/${document.location.pathname.split('/')[1]}_library.xlsx`, {
          method: 'HEAD' 
        });
        if (response.ok) { 
          setFileExists(true);
        }
      } catch (error) {
        console.error('Error checking file existence:', error);
      }
    };

    checkFileExists();
  }, []);

  return (
    <Grid item>
      {fileExists && 
        <IconLink IconElement={DownloadIcon} text="Download" isDownload={true} link={`${document.location.pathname.split('/')[1]}/data/${document.location.pathname.split('/')[1]}_library.xlsx`}></IconLink>}
    </Grid>
  );
};

export default CheckFileAndDownloadIcon;
