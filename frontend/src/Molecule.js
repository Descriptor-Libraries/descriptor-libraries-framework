import React from "react";
import CircularProgress from "@mui/material/CircularProgress";
import Typography from "@mui/material/Typography";
import { makeStyles } from "@mui/material/styles";
import ButtonGroup from "@mui/material/ButtonGroup";
import Button from "@mui/material/Button";
import Grid from "@mui/material/Grid";
import { useParams, useHistory, useLocation } from "react-router-dom";

import PubChemInfo from "./PubChem.js";
import CDKDepict from "./CDKDepict.js";
import MoleculeDataGrid from "./MoleculeDataGrid.js";

const useStyles = makeStyles({
  root: {
    flexGrow: 1,
  },
});

function Molecule(props) {
  const [data, setData] = React.useState(null);
  const [error, setError] = React.useState(null);
  const [loading, setLoading] = React.useState(true);
  let { moleculeId } = useParams();

  const fetchData = () => {
    setError(null);
    fetch("/api/v1/molecule/" + encodeURIComponent(moleculeId))
      .then((r) => {
        if (r.status == 200) return r.json();
        setError(`${r.status} ${r.statusText}`);
        return {};
      })
      .then((data) => {
        setData(data);
      })
      .finally(() => setLoading(false));
  };

  React.useEffect(fetchData, [moleculeId]);

  if (loading) return <CircularProgress />;
  if (error) {
    return (
      <div>
        <Typography variant="h6">
          Oh no! Something went wrong with Molecule #{moleculeId}! :'(
        </Typography>
        <Typography>{error}</Typography>
      </div>
    );
  }
  return (
    <div>
      <Grid
        container
        direction="row"
        justify="space-between"
        alignItems="center"
      >
        <Grid item>
          <Typography variant="h4">Molecule #{moleculeId}</Typography>
        </Grid>
        <Grid item>
          <Button
            variant="contained"
            color="primary"
            href={`data:text/json;charset=utf-8,${encodeURIComponent(
              JSON.stringify(data)
            )}`}
            download={`molecule-${moleculeId}.json`}
          >
            Download
          </Button>
        </Grid>
      </Grid>
      <MoleculeData data={data} />
    </div>
  );
}

function MoleculeData(props) {
  let history = useHistory();
  let currentLocation = useLocation();

  const handleClickConformer = () => {
    history.push("/conformer/" + props.data.conformers_id[0], {
      prevPath: currentLocation,
    });
  };
  return (
    <div>
      <CDKDepict width={300} height={300} zoom={8} smiles={props.data.smiles} />
      <PubChemInfo smiles={props.data.smiles} />
      <Typography variant="h6">Conformers</Typography>
      {props.data.dft_data != undefined ? (
        <ButtonGroup variant="text">
          <Button onClick={handleClickConformer}>
            See conformer structures
          </Button>
        </ButtonGroup>
      ) : (
        <Typography variant="body">
          No conformers have been computed for this molecule
        </Typography>
      )}
      <MoleculeDataGrid data={props.data} />
    </div>
  );
}

export default Molecule;
