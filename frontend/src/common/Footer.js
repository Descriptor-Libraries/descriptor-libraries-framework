import React from 'react';
import Box from '@mui/material/Box';
import Grid from '@mui/material/Grid';
import Divider from '@mui/material/Divider';

import ccas_logo from "../images/CCAS_logo_light.png";
import molssi_logo from "../images/molssi_main_logo.png";

import "../assets/custom.css";

function Footer() {
    return (
    <><Box class="footer">
        <Divider />
        <Grid container spacing={2}>
            <Grid item xs={2}>
                <a href="https://ccas.nd.edu/" target="_blank" title="Go to C-CAS in a new tab" >
                <img src={ccas_logo} alt="CCAS logo" class="footer_logo" />
                </a>
            </Grid>
            <Grid item xs={8}>
                <p style={{ textAlign: 'left' }}> &copy; Copyright 2019-2023 <a href="https://molssi.org/">The Molecular Sciences Software Institute </a>
                    and the <a href="https://ccas.nd.edu/">Center for Computer Assisted Synthesis</a>
                </p>
                <p style={{ textAlign: 'left' }}>
                    MolSSI is Funded by the National Science Foundation <a href="https://nsf.gov/awardsearch/showAward?AWD_ID=1547580">OAC-1547580 </a> and
                    <a href="https://www.nsf.gov/awardsearch/showAward?AWD_ID=2136142"> CHE-2136142.</a>
                </p>
            </Grid>
            <Grid item xs={2}>
                <a href="https://molssi.org/" target="_blank" title="Go to MolSSI in a new tab">
                    <img src={molssi_logo} alt="Go to molssi.org" class="footer_logo"/>
                </a>
            </Grid>
        </Grid>
    </Box>
    </>
    )
    }

    export default Footer;