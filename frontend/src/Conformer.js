import React from "react";
import ArrowBackIosIcon from "@material-ui/icons/ArrowBackIos";
import Button from "@material-ui/core/Button";
import CircularProgress from "@material-ui/core/CircularProgress";
import Grid from "@material-ui/core/Grid";
import Typography from "@material-ui/core/Typography";
import { useParams, useHistory } from "react-router-dom";
import * as ngl from "ngl";

function Conformer(props) {
  const [error, setError] = React.useState(null);
  const [othersId, setOthersId] = React.useState([]);
  let { conformerId } = useParams();
  let history = useHistory();
  const handleBackButton = () => {
    console.log(history.location.state);
    if (
      history.location.state != undefined &&
      history.location.state.prevPath != undefined
    ) {
      history.push(history.location.state.prevPath.pathname);
    } else {
      history.goBack();
    }
  };

  const handleClickConformer = (confId) => {
    history.push("/conformer/" + confId, history.location.state);
  };

  React.useEffect(() => {
    fetch("/api/v1/conformer/others_id/" + conformerId)
      .then((res) => {
        if (res.status == 200) return res.json();
        setError(res);
        return [];
      })
      .then((data) => setOthersId(data));
  }, []);
  return (
    <Grid container direction="row" justify="center" alignItems="flex-start">
      <Grid item xs={12} sm={6} xl={12}>
        <Button startIcon={<ArrowBackIosIcon />} onClick={handleBackButton}>
          Back
        </Button>
      </Grid>
      <Grid item xs={12} sm={8} xl={6}>
        <iframe
          src={"/nglviewer.html?conformerId=" + conformerId}
          width="850px"
          height="850px"
          style={{ border: 0 }}
        />
      </Grid>
      <Grid item xs={12} sm={8} xl={6}>
        <Typography variant="h6">
          Other conformers of the same molecule
        </Typography>
        {othersId.map((confId) => (
          <Button onClick={() => handleClickConformer(confId)}>{confId}</Button>
        ))}
      </Grid>
    </Grid>
  );
}

export default Conformer;
