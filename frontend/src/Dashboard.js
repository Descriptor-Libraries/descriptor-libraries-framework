import React, { lazy, Suspense } from "react";
import CircularProgress from "@mui/material/CircularProgress";
import Button from "@mui/material/Button";
import Popover from "@mui/material/Popover";
import Grid from "@mui/material/Grid";
import FormGroup from "@mui/material/FormGroup";
import FormControlLabel from "@mui/material/FormControlLabel";
import FormSwitch from "@mui/material/Switch";
import FormLabel from "@mui/material/FormLabel";
import RadioGroup from "@mui/material/RadioGroup";
import Radio from "@mui/material/Radio";
import Tour from "reactour";
import Typography from "@mui/material/Typography";
import { makeStyles } from "@mui/material/styles";
import * as dfd from "danfojs/src/index";

import {
  useParams,
  useRouteMatch,
  useHistory,
  Switch,
  Route,
} from "react-router-dom";

import Plot from "react-plotly.js";
import DrawSmiles from "./DrawSmiles.js";

const Molecule = lazy(() => import("./Molecule.js"));

const useStyles = makeStyles((theme) => ({
  plot: {
    marginTop: theme.spacing(2),
  },
  popover: {
    pointerEvents: "none",
  },
  paper: {},
}));

function PCAPlot(props) {
  const classes = useStyles();

  const [data, setData] = React.useState([]);
  const [dataframe, setDataframe] = React.useState(null);

  const [showMlData, setShowMlData] = React.useState(false);

  const [visibleGroups, setVisibleGroups] = React.useState([
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
    true,
  ]);
  const [proj, setProj] = React.useState("umap");

  const labels = [
    "P[C]<sub>3</sub>",
    "P[C]<sub>n</sub>[R]<sub>m</sub>",
    "P[O]<sub>3</sub>",
    "P[O]<sub>x</sub>[N]<sub>y</sub>",
    "P[N]<sub>3</sub>",
    "P[C]<sub>n</sub>[N]<sub>m</sub>",
    "PF<sub>n</sub>[R]<sub>m</sub>",
    "P[Si]<sub>n</sub>[R]<sub>m</sub>",
    "Others",
  ];

  const palette = [
    "#ce3042",
    "#e3e300",
    "#0067a5",
    "#f38400",
    "#a756b2",
    "#8db600",
    "#a1caf1",
    "#848482",
    "#654522",
  ];

  const toggleVisibleGroup = (group) => {
    let newVisibleGroups = visibleGroups.map((e, i) => {
      if (group == i) return !e;
      return e;
    });
    setVisibleGroups(newVisibleGroups);
  };
  const fetchData = () => {
    dfd.read_csv("https://kraken.cs.toronto.edu/ligand_data.csv").then((df) => {
      setDataframe(df);
    });
  };

  const df2plotly = (df, color, name) => {
    return {
      x: df[`${proj}1`].data,
      y: df[`${proj}2`].data,
      name: name,
      mode: "markers",
      marker: { opacity: 0.8, color: color, symbol: "square", size: 3 },
      type: "scatter",
      customdata: df.data.map((e) => ({ molecule_id: e[0], smiles: e[1] })),
    };
  };

  const formatData = () => {
    let df = dataframe;
    if (dataframe === null) {
      setData([]);
      return;
    }
    if (showMlData === false) {
      df = dataframe.query({ column: "ID", is: "<", to: 1545 });
    }
    let plotlyData = visibleGroups.map((e, i) => {
      if (e) {
        let op = "==";
        if (i > 7) op = ">=";
        return df2plotly(
          df.query({ column: "pat_enc", is: op, to: i }),
          palette[i],
          labels[i]
        );
      }
      return false;
    });
    setData(
      plotlyData.filter((e) => {
        return e != false;
      })
    );
  };

  const [layout, setLayout] = React.useState({
    hovermode: "closest",
    width: window.innerHeight,
    height: window.innerHeight * 0.8,
    xaxis: {
      title: { text: `${proj} 1` },
    },
    yaxis: {
      title: { text: `${proj} 2` },
    },
    legend: {
      y: 0.5,
      yref: "paper",
    },
  });

  React.useEffect(fetchData, []);
  React.useEffect(formatData, [dataframe, visibleGroups, proj, showMlData]);

  return (
    <div>
      <Plot
        className={classes.plot}
        data={data}
        layout={layout}
        onUpdate={(figure) => {
          setLayout(figure.layout);
        }}
        onClick={(e) => {
          props.setMoleculeId(e.points[0].customdata.molecule_id);
        }}
        onHover={(e) => {
          props.openPopover(
            e.points[0].customdata.smiles,
            e.event.y - 10,
            e.event.x - 10
          );
        }}
        onUnhover={(e) => {
          props.closePopover(false);
        }}
      />
      <div>
        <FormLabel component="legend">Data Source</FormLabel>
        <FormGroup row>
          <FormControlLabel
            control={
              <FormSwitch
                checked={showMlData}
                onChange={() => setShowMlData(!showMlData)}
                name="Show ML data"
                color="primary"
              />
            }
            label={<Typography>Show ML Data</Typography>}
          />
        </FormGroup>
      </div>
      <div data-tut="projection">
        <FormLabel component="legend">Projection</FormLabel>
        <FormGroup row>
          <RadioGroup
            row
            aria-label="projection"
            name="projection"
            value={proj}
            onChange={(e) => {
              setProj(e.target.value);
            }}
          >
            <FormControlLabel value="pca" control={<Radio />} label="PCA" />
            <FormControlLabel value="umap" control={<Radio />} label="UMAP" />
          </RadioGroup>
        </FormGroup>
      </div>
      <div data-tut="groups">
        <FormLabel component="legend">Data Visibility</FormLabel>
        <FormGroup row>
          {visibleGroups.map((e, i) => (
            <FormControlLabel
              control={
                <FormSwitch
                  checked={visibleGroups[i]}
                  onChange={() => toggleVisibleGroup(i)}
                  name={labels[i]}
                  color="primary"
                />
              }
              label={
                <Typography dangerouslySetInnerHTML={{ __html: labels[i] }} />
              }
            />
          ))}
        </FormGroup>
      </div>
    </div>
  );
}

function HoverPoint(props) {
  const classes = useStyles();

  return (
    <Popover
      disableScrollLock={true}
      className={classes.popover}
      open={props.open}
      anchorReference="anchorPosition"
      anchorPosition={{ top: props.anchorTop, left: props.anchorLeft }}
      anchorOrigin={{
        vertical: "bottom",
        horizontal: "right",
      }}
      transformOrigin={{
        vertical: "bottom",
        horizontal: "right",
      }}
      onClose={() => {
        props.closePopover(false);
      }}
      elevation={3}
      disableRestoreFocus
    >
      <DrawSmiles smiles={props.smiles} width={150} height={150} />
    </Popover>
  );
}

const tourSteps = [
  {
    content:
      "Welcome to the Kraken WebApp! Kraken stands for Kolossal viRtual dAtabase moleKular dEscriptors of orgaNophosphorus ligands. This WebApp give you access to a database of monodentate organophosphorus(iii) ligands.",
  },
  {
    content:
      "To understand the data we generated better, we computed descriptors of the DFT-computed compounds. This plot is a 2D projection of the descriptors we generated.",
    selector: '[data-tut="pca__plot"]',
  },
  {
    content:
      "You can chose among 2 different projection: principal component analysis (PCA) and uniform manifold approximation and projection (UMAP)",
    selector: '[data-tut="projection"]',
  },
  {
    content:
      "You can also chose to enable/disable different group of ligands, according to their structures",
    selector: '[data-tut="groups"]',
  },
  {
    content:
      "If you click on a point of the graph, further information are displayed here",
    selector: '[data-tut="molecule"]',
  },
];

function Dashboard() {
  let { moleculeId } = useParams();
  let match = useRouteMatch();
  let history = useHistory();
  const [smiles, setSmiles] = React.useState("");
  const [anchorTop, setAnchorTop] = React.useState(0);
  const [anchorLeft, setAnchorLeft] = React.useState(0);
  const [open, setOpen] = React.useState(false);
  const [isTourOpen, setIsTourOpen] = React.useState(false);

  const openPopover = (s, t, l) => {
    setSmiles(s);
    setAnchorTop(t);
    setAnchorLeft(l);
    setOpen(true);
  };

  const handlePointClicked = (molId) => {
    history.push(`/dashboard/${molId}`);
  };

  return (
    <div>
      <Tour
        steps={tourSteps}
        isOpen={isTourOpen}
        onRequestClose={() => setIsTourOpen(false)}
        accentColor="#192b28"
        rounded={5}
      />
      <Grid container direction="column" justify="center" alignItems="center">
        <Grid item>
          <Typography variant="h3">Projection of the Chemical Space</Typography>
          <Typography variant="body1">
            To keep the visualisation fluid, only 4545 data points (ie
            molecules) have been loaded. There is in total 330949 data points
            available in the database.
          </Typography>
          <Typography variant="body1">
            Not sure how to use this page?
          </Typography>
          <Button
            variant="contained"
            color="primary"
            onClick={() => {
              history.push("/dashboard/1");
              setIsTourOpen(true);
            }}
            disableElevation
          >
            Take the tour!
          </Button>
        </Grid>
        <Grid item>
          <div data-tut="pca__plot">
            <Grid
              container
              direction="row"
              justify="center"
              alignItems="stretch"
            >
              <PCAPlot
                setMoleculeId={handlePointClicked}
                openPopover={openPopover}
                closePopover={setOpen}
              />
            </Grid>
          </div>
        </Grid>
        <Grid item xs={12} sm={12} md={8} xl={6}>
          <Suspense fallback={<CircularProgress />}>
            <Switch>
              <Route path={`${match.path}/:moleculeId`}>
                <div data-tut="molecule">
                  <Molecule moleculeId={moleculeId} />
                </div>
              </Route>
              <Route>
                <Typography variant="h6">
                  Click on a point of the plot!
                </Typography>
              </Route>
            </Switch>
          </Suspense>
        </Grid>
        <HoverPoint
          open={open}
          smiles={smiles}
          anchorTop={anchorTop}
          anchorLeft={anchorLeft}
        />
      </Grid>
    </div>
  );
}

export default Dashboard;
