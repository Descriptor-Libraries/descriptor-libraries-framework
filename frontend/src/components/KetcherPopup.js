import React, { useEffect, useState } from 'react';

import Button from '@mui/material/Button';
import Dialog from '@mui/material/Dialog';
import Slide from '@mui/material/Slide';
import Box from '@mui/material/Box';

import { Ketcher } from 'ketcher-core';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Editor } from 'ketcher-react';
import "ketcher-react/dist/index.css";

const structServiceProvider = new StandaloneStructServiceProvider();

const Transition = React.forwardRef(function Transition(props, ref) {
  return <Slide direction="up" ref={ref} {...props} />;
});

export default function FullScreenDialog({ representation, smilesChange, smartsChange }) {
  const [ open, setOpen ] = React.useState(false);
  const [ smiles, setSmiles ] = useState();
  const [ SMARTS, setSMARTS ] = useState();
  const [ ketcher, setKetcher ] = useState();
    
  const handleClickOpen = () => {
      setOpen(true);
  };

  const handleClose = () => {
      if (representation === "smiles"){
          ketcher.getSmiles().then(result => {setSmiles(result);});
      }
      if (representation === "SMARTS"){
          ketcher.getSmarts().then(result => {setSMARTS(result);});
      }
      setOpen(false);
  };

  useEffect(() => {
    smartsChange(SMARTS);
  }, [SMARTS]);

  useEffect(() => {
    smilesChange(smiles);
  }, [smiles]);

  return (
    <div>
      <Button onClick={handleClickOpen}>
        DRAW
      </Button>
      <Dialog
        fullWidth={true}
        maxWidth={"lg"}
        open={open}
        onClose={handleClose}
        TransitionComponent={Transition}
      >
          <Editor
              staticResourcesUrl={process.env.PUBLIC_URL}
              structServiceProvider={structServiceProvider}
              onInit={(ketcher) => {
                  setKetcher(ketcher)
              }}
          />
      </Dialog>
    </div>
  );
}