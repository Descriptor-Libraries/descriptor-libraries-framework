import React from "react";
import Typography from "@material-ui/core/Typography";
import Button from "@material-ui/core/Button";
import Table from "@material-ui/core/Table";
import TableBody from "@material-ui/core/TableBody";
import TableCell from "@material-ui/core/TableCell";
import TableContainer from "@material-ui/core/TableContainer";
import TableRow from "@material-ui/core/TableRow";
import Paper from "@material-ui/core/Paper";
import CircularProgress from "@material-ui/core/CircularProgress";

export default function PubChemInfo(props) {
  const baseURL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/";
  let [error, setError] = React.useState(null);
  let [loading, setLoading] = React.useState(true);
  let [iupacName, setIupacName] = React.useState("");
  let [inchi, setInchi] = React.useState("");
  let [cid, setCid] = React.useState(-1);

  React.useEffect(() => {
    setLoading(true);
    fetch(
      baseURL +
        "compound/smiles/" +
        encodeURIComponent(props.smiles) +
        "/record/json"
    )
      .then((r) => {
        if (r.status == 200) {
          return r.json();
        }
        setError(r);
        return {};
      })
      .then((data) => {
        if (data.hasOwnProperty("PC_Compounds") == false) {
          setError("Compound not found on PubChem!");
          return;
        }
        const props = data.PC_Compounds[0].props;
        if (props[6] != undefined) setIupacName(props[6].value.sval);
        else setIupacName("N/A");
        if (props[12] != undefined) setInchi(props[12].value.sval);
        else setInchi("N/A");
        if (data.PC_Compounds[0].id.id != undefined)
          setCid(data.PC_Compounds[0].id.id.cid);
        else setCid("N/A");
      })
      .finally(() => {
        setLoading(false);
      });
  }, [props.smiles]);

  if (loading)
    return (
      <div>
        <Typography variant="h6">PubChem Data</Typography>
        <CircularProgress />
      </div>
    );
  if (error)
    return (
      <div>
        <Typography variant="h6">PubChem Data</Typography>
        <p>Sorry, fetching data on PubChem did not work ¯\_(ツ)_/¯: {error}</p>
      </div>
    );

  return (
    <div>
      <Typography variant="h6">PubChem Data</Typography>
      <TableContainer component={Paper}>
        <TableBody>
          <TableRow key="IUPAC Name">
            <TableCell align="right">IUPAC Name</TableCell>
            <TableCell aligh="right">{iupacName}</TableCell>
          </TableRow>
          <TableRow key="InchI">
            <TableCell align="right">InchI</TableCell>
            <TableCell aligh="right">{inchi}</TableCell>
          </TableRow>
          <TableRow key="CID">
            <TableCell align="right">CID</TableCell>
            <TableCell aligh="right">{cid}</TableCell>
          </TableRow>
        </TableBody>
      </TableContainer>
      <Button
        href={
          "https://pubchem.ncbi.nlm.nih.gov/#query=" +
          props.smiles +
          "&input_type=smiles"
        }
      >
        Lookup on PubChem
      </Button>
    </div>
  );
}
