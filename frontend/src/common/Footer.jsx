import React, { useEffect, useState } from 'react';
import Box from '@mui/material/Box';
import Grid from '@mui/material/Grid';
import Divider from '@mui/material/Divider';

import ccas_logo from "../images/CCAS_logo_light.png";
import molssi_logo from "../images/molssi_main_logo.png";

import "../assets/custom.css";

async function fetchVersion(url) {
    /**
     * Reads a version document: {"version": ..., "git_sha": ..., "build_date": ...}.
     * @param {string} url Where to read it from.
     * @return {object|null} The document, or null when it is not available
     * (for example the dev server has no version.json, and older backends have no /version).
     */
    try {
        const response = await fetch(url);
        if (!response.ok) {
            return null;
        }
        const data = await response.json();
        return data && data.version ? data : null;
    } catch (error) {
        return null;
    }
}

function versionLabel(info) {
    /**
     * Formats a version document for display, e.g. "v2026.11.0 (c6e4675)".
     * @param {object|null|undefined} info The version document; undefined while it is loading.
     * @return {string} The text to show.
     */
    if (info === undefined) {
        return "…";
    }
    if (!info) {
        return "unknown";
    }
    const version = /^\d/.test(info.version) ? `v${info.version}` : info.version;
    const sha = info.git_sha && info.git_sha !== "unknown" ? ` (${info.git_sha.slice(0, 7)})` : "";
    return `${version}${sha}`;
}

function Footer() {
    const [isMobile, setIsMobile] = useState(window.innerWidth < 768);
    const [frontendVersion, setFrontendVersion] = useState(undefined);
    const [backendVersion, setBackendVersion] = useState(undefined);

    useEffect(() => {
        // The library is the first part of the path, as everywhere else in the app.
        const library = document.location.pathname.split('/')[1];
        // version.json is written into the frontend image when it is built.
        fetchVersion(`/${library}/version.json`).then(setFrontendVersion);
        fetchVersion(`/api/${library}/version`).then(setBackendVersion);
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


    return (
    <><Box className="footer">
        <Divider />
        <Grid container spacing={2} alignItems="center" justifyContent="center">
            
        <Grid item xs={2} style={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                {!isMobile &&
                    <a href="https://ccas.nd.edu/" target="_blank" rel="noreferrer" title="Go to C-CAS in a new tab" 
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
                <p className="footer_version" style={{ textAlign: 'center' }}
                    title={[frontendVersion && `Frontend built ${frontendVersion.build_date}`, backendVersion && `Backend built ${backendVersion.build_date}`].filter(Boolean).join(' · ')}>
                    Frontend {versionLabel(frontendVersion)} · Backend {versionLabel(backendVersion)}
                </p>

            </Grid>
            <Grid item xs={2}>
            {!isMobile &&
                <a href="https://molssi.org/" target="_blank" rel="noreferrer" title="Go to MolSSI in a new tab">
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