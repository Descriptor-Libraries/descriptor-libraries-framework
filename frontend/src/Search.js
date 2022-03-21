import React from "react";
import { useTheme } from "@material-ui/core/styles";
import useMediaQuery from "@material-ui/core/useMediaQuery";
import Button from "@material-ui/core/Button";
import Card from "@material-ui/core/Card";
import CardContent from "@material-ui/core/CardContent";
import CardActions from "@material-ui/core/CardActions";
import CircularProgress from "@material-ui/core/CircularProgress";
import Grid from "@material-ui/core/Grid";
import FormControlLabel from "@material-ui/core/FormControlLabel";
import FormGroup from "@material-ui/core/FormGroup";
import Switch from "@material-ui/core/Switch";
import Typography from "@material-ui/core/Typography";
import CDKDepict from "./CDKDepict.js";
import { useHistory } from "react-router-dom";
import Tour from "reactour";
import KetcherDialog from "./KetcherDialog.js";

function SearchHeader(props) {
  const [openKetcher, setOpenKetcher] = React.useState(false);
  return (
    <div>
      <Grid
        container
        style={{ marginTop: "10px" }}
        classdirection="column"
        justify="space-between"
        alignItems="center"
      >
        <Grid item>
          <div data-tut="structure">
            <Typography variant="h6">Looking for:</Typography>
            <CDKDepict
              smiles={props.smiles}
              width={150}
              height={150}
              zoom={8}
            />
          </div>
          <div data-tut="edit">
            <Button onClick={() => setOpenKetcher(true)}>Edit</Button>
          </div>
        </Grid>
        <div data-tut="data">
          <Grid item xs={12} sm={12} md={6} lg={6} xl={4}>
            <FormGroup row>
              <FormControlLabel
                control={
                  <Switch
                    name="dft"
                    checked={props.switchState.dft}
                    onChange={props.handleChange}
                  />
                }
                label="DFT"
              />
              <FormControlLabel
                control={
                  <Switch
                    name="xtb"
                    checked={props.switchState.xtb}
                    onChange={props.handleChange}
                  />
                }
                label="XTB"
              />
              <FormControlLabel
                control={
                  <Switch
                    name="ml"
                    checked={props.switchState.ml}
                    onChange={props.handleChange}
                  />
                }
                label="ML"
              />
            </FormGroup>
          </Grid>
        </div>
      </Grid>
      <KetcherDialog
        open={openKetcher}
        onClose={(smiles) => {
          props.setSmiles(smiles);
          setOpenKetcher(false);
        }}
      />
    </div>
  );
}

function SearchResultCard(props) {
  let history = useHistory();
  let theme = useTheme();

  const sm = useMediaQuery(theme.breakpoints.up("sm"));

  return (
    <Card elevation={3}>
      <CardContent>
        <CDKDepict
          smiles={props.data.smiles}
          width={sm ? 150 : 250}
          height={sm ? 100 : 200}
          zoom={sm ? 2.5 : 4}
        />
      </CardContent>
      <CardActions>
        <Button
          size="small"
          onClick={() => history.push("/molecule/" + props.data.molecule_id)}
        >
          View
        </Button>
        <Button
          size="small"
          onClick={() => navigator.clipboard.writeText(props.data.smiles)}
        >
          Copy SMILES
        </Button>
      </CardActions>
    </Card>
  );
}

function SearchResults(props) {
  if (props.loading) {
    return (
      <Grid item xs={10} sm={10} md={8} xl={8}>
        <CircularProgress />
      </Grid>
    );
  }
  if (props.error) {
    console.log(props.error);
    return <Typography variant="h6">Oh no! Something went wrong!</Typography>;
  }

  return (
    <Grid
      container
      spacing={2}
      direction="row"
      justify="center"
      alignItems="flex-start"
    >
      {props.data.map((elem) => (
        <Grid item>
          <SearchResultCard data={elem} />
        </Grid>
      ))}
    </Grid>
  );
}

function Search(props) {
  const [smiles, setSmiles] = React.useState("P");
  const [src, setSrc] = React.useState({
    xtb: true,
    dft: true,
    ml: false,
  });
  const [data, setData] = React.useState([]);
  const [loading, setLoading] = React.useState(true);
  const [error, setError] = React.useState(false);
  const [isTourOpen, setIsTourOpen] = React.useState(false);

  const changeSrc = (event) => {
    setSrc({ ...src, [event.target.name]: event.target.checked });
  };

  const steps = [
    {
      content: "Here is the current strucutre we are looking for",
      selector: '[data-tut="structure"]',
    },
    {
      content: "You can change it using this button that will open an editor.",
      selector: '[data-tut="edit"]',
    },
    {
      content: "You can also select the method used to compute the data.",
      selector: '[data-tut="data"]',
    },
  ];

  const loadData = (limit, append) => {
    let skip = data.length == undefined ? 0 : data.length;
    fetch(
      "/api/v1/molecule/search" +
        "?substructure=" +
        encodeURIComponent(smiles) +
        "&with_dft=" +
        src.dft +
        "&with_xtb=" +
        src.xtb +
        "&with_ml=" +
        src.ml +
        "&skip=" +
        skip +
        "&limit=" +
        limit,
      { method: "PUT", mode: "cors" }
    )
      .then((res) => {
        if (res.status == 200) return res.json();
        setError(res);
        return {};
      })
      .then((d) => {
        if (append) setData([...data, ...d]);
        else setData(d);
      })
      .finally(setLoading(false));
  };
  React.useEffect(() => loadData(20, false), [smiles]);
  return (
    <div>
      <Tour
        steps={steps}
        isOpen={isTourOpen}
        onRequestClose={() => setIsTourOpen(false)}
        accentColor="#192b28"
        rounded={5}
      />
      <Grid container direction="column" justify="center" alignItems="center">
        <Grid item>
          <Typography variant="h3">Substructure Search</Typography>
          <Typography variant="body1">
            Not sure how to use this page?
          </Typography>
          <Button
            variant="contained"
            color="primary"
            onClick={() => {
              setIsTourOpen(true);
            }}
            disableElevation
          >
            Take the tour!
          </Button>
        </Grid>
        <Grid item>
          <Grid container direction="row" justify="center" alignItems="stretch">
            <Grid item xs={10} sm={10} md={8} xl={6}>
              <SearchHeader
                switchState={src}
                handleChange={changeSrc}
                smiles={smiles}
                setSmiles={setSmiles}
              />
            </Grid>
          </Grid>
          <Grid item xs={12} sm={12} md={12} xl={12}>
            <SearchResults data={data} loading={loading} error={error} />
          </Grid>
          <Grid container justify="center" item xs={10} sm={12} md={10} xl={10}>
            <Button
              variant="contained"
              color="secondary"
              style={{ margin: "10px" }}
              onClick={() => {
                loadData(20, true);
              }}
            >
              Load More
            </Button>
          </Grid>
        </Grid>
      </Grid>
    </div>
  );
}

export default Search;
