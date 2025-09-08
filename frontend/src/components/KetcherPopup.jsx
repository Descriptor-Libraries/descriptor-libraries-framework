import React, { useState, useEffect } from 'react';

import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import Slide from '@mui/material/Slide';
import CloseIcon from '@mui/icons-material/Close';
import CheckIcon from '@mui/icons-material/Check';
import IconButton from '@mui/material/IconButton';
import AppBar from '@mui/material/AppBar';
import Toolbar from '@mui/material/Toolbar';
import Typography from '@mui/material/Typography';

import { ChemicalMimeType } from 'ketcher-core';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Editor } from 'ketcher-react';
import "ketcher-react/dist/index.css";

// Async initialization function
async function getStructServiceProvider() {
  return new StandaloneStructServiceProvider();
}

const Transition = React.forwardRef(function Transition(props, ref) {
  return <Slide direction="up" ref={ref} {...props} />;
});

export default function FullScreenDialog({ ketcherCallBack }) {
  const [open, setOpen] = React.useState(false);
  const [ketcher, setKetcher] = useState();
  const [isReady, setIsReady] = useState(false);
  const [applying, setApplying] = useState(false);
  
  // Proper async initialization
  const [structServiceProvider, setStructServiceProvider] = useState(null);
  
  useEffect(() => {
    getStructServiceProvider().then(setStructServiceProvider);
  }, []);

  const handleClickOpen = () => setOpen(true);
  const handleCancel = () => setOpen(false);

  async function handleApply() {
    if (!ketcher || applying) return;
    setApplying(true);
    try {
      let smiles = '';
      let smarts = '';
      // Prefer direct SMILES; normalize with Indigo if available
      try {
        smiles = await ketcher.getSmiles();
      } catch (e) {
        // no-op, try Indigo path below
      }
      try {
        if (ketcher.indigo) {
          const mol = await ketcher.indigo.aromatize(smiles || (await ketcher.getMolfile()));
          const converted = await ketcher.indigo.convert(mol, { outputFormat: ChemicalMimeType.DaylightSmiles });
          smiles = converted?.struct || smiles;
        }
      } catch (e) {
        // keep whatever we have in smiles
      }
      try {
        smarts = await ketcher.getSmarts();
      } catch (e) {
        // best-effort only
      }
      ketcherCallBack([smiles, smarts]);
    } finally {
      setApplying(false);
      setOpen(false);
    }
  }

  // Show loading state while service provider initializes
  if (!structServiceProvider) {
    return (
      <div>
        <Button variant="contained" sx={{ mb: 5 }} disabled>
          <span style={{ textTransform: 'capitalize', fontSize: '16px' }}>
            Loading Molecular Sketcher...
          </span>
        </Button>
      </div>
    );
  }

  return (
    <div>
      <Button variant="contained" sx={{ mb: 5 }} onClick={handleClickOpen}>
        <span style={{ textTransform: 'capitalize', fontSize: '16px' }}>
          Open Molecular Sketcher
        </span>
      </Button>
      <Dialog
        fullScreen
        maxWidth={"lg"}
        open={open}
        onClose={handleCancel}
        keepMounted
        TransitionComponent={Transition}
      >
        <AppBar sx={{ position: 'relative' }}>
          <Toolbar>
            <IconButton
              edge="start"
              color="inherit"
              onClick={handleCancel}
              aria-label="close"
            >
              <CloseIcon />
            </IconButton>
            <Button
              color="inherit"
              startIcon={<CheckIcon />}
              onClick={handleApply}
              disabled={!isReady || applying}
              sx={{ ml: 1 }}
            >
              Use structure
            </Button>
          </Toolbar>
        </AppBar>
        <Editor
          staticResourcesUrl={process.env.PUBLIC_URL}
          structServiceProvider={structServiceProvider}
          onInit={(inst) => {
            window.ketcher = inst;
            setKetcher(inst);
            setIsReady(true);
            // optional: ensure canvas sizes correctly on first paint
            setTimeout(() => { try { inst.editor?.resize?.(); } catch {} }, 0);
          }}
        />
      </Dialog>
    </div>
  );
}