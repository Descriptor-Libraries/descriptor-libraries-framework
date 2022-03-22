import Button from "@mui/material/Button";
import Dialog from "@mui/material/Dialog";
import DialogTitle from "@mui/material/DialogTitle";
import DialogContent from "@mui/material/DialogContent";
import DialogActions from "@mui/material/DialogActions";

import { Editor } from "ketcher-react";
import { StandaloneStructServiceProvider } from "ketcher-standalone";
import "ketcher-react/dist/index.css";

function KetcherDialog(props) {
  const handleClick = () => {
    window.ketcher.getSmilesAsync().then((smiles) => {
      if (smiles == "") props.onClose("P");
      else props.onClose(smiles);
    });
  };
  let structServiceProvider = new StandaloneStructServiceProvider();
  return (
    <Dialog fullWidth maxWidth={"xl"} open={props.open} onClose={handleClick}>
      <DialogTitle>Edit Structure</DialogTitle>
      <DialogContent>
        <div style={{ height: 600 }}>
          <Editor
            style={{ height: 600 }}
            structServiceProvider={structServiceProvider}
          />
        </div>
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClick}>Done</Button>
      </DialogActions>
    </Dialog>
  );
}

export default KetcherDialog;
