import React, { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import Grid from '@mui/material/Grid';
import Divider from '@mui/material/Divider';

import ccas_logo from "../images/CCAS_logo_light.png";
import molssi_logo from "../images/molssi_main_logo.png";

import "../assets/custom.css";

function Footer() {
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


    return (
    <><Box className="footer">
        <Divider />
        <Grid container spacing={2} alignItems="center" justifyContent="center">
            
        <Grid item xs={2} style={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                {!isMobile &&
                    <a href="https://ccas.nd.edu/" target="_blank" title="Go to C-CAS in a new tab" 
                    style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', width: '100%', height: '100%' }}>
                        <img src={ccas_logo} alt="CCAS logo" className="footer_logo" />
                    </a>
                }
            </Grid>
            <Grid item xs={8}>
                <p style={{ textAlign: 'center' }}>This site is maintained by <a href="https://molssi.org/">The Molecular Sciences Software Institute (MolSSI)</a>
                    and the <a href="https://ccas.nd.edu/">Center for Computer Assisted Synthesis (CCAS).</a>
                </p>
                <p style={{ textAlign: 'center' }}>
                    MolSSI is Funded by the National Science Foundation <a href="https://nsf.gov/awardsearch/showAward?AWD_ID=1547580">OAC-1547580 </a> and
                    <a href="https://www.nsf.gov/awardsearch/showAward?AWD_ID=2136142"> CHE-2136142.</a> CCAS is Funded by the National Science Foundation <a href="https://www.nsf.gov/awardsearch/showAward?AWD_ID=2202693">CHE–2202693.</a>
                </p>
                
            </Grid>
            <Grid item xs={2}>
            {!isMobile &&
                <a href="https://molssi.org/" target="_blank" title="Go to MolSSI in a new tab">
                    <img src={molssi_logo} alt="Go to molssi.org" className="footer_logo"/>
                </a>
             }
            </Grid>

        </Grid>
    </Box>
    </>
    )
    }

    export default Footer;