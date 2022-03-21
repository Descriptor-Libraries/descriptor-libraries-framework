import {
  BrowserRouter as Router,
  Switch,
  Route,
  useHistory,
} from "react-router-dom";
import { Suspense, lazy } from "react";
import ArrowBackIosIcon from "@mui/icons-material/ArrowBackIos";
import CircularProgress from "@mui/material/CircularProgress";
import Button from "@mui/material/Button";
import Grid from "@mui/material/Grid";
import { createMuiTheme, ThemeProvider } from "@mui/material/styles";

import About from "./About.js";
import Contact from "./Contact.js";
import Download from "./Download.js";
import TopBar from "./TopBar.js";
import "./App.css";

const Search = lazy(() => import("./Search.js"));
const Conformer = lazy(() => import("./Conformer.js"));
const Molecule = lazy(() => import("./Molecule.js"));
const Dashboard = lazy(() => import("./Dashboard.js"));

const theme = createMuiTheme({
  typography: {
    h6: {
      fontFamily: "Helvetica",
    },
  },
  palette: {
    primary: {
      contrastText: "#12f98b",
      main: "#192b28",
    },
    secondary: {
      main: "#ffeb00",
    },
    mtYellow: {
      main: "#ffeb00",
    },
    mtPink: {
      main: "#e30060",
    },
    mtBlue: {
      main: "#0095e3",
    },
    background: {
      paper: "#fff",
    },
  },
});

function MoleculePage() {
  let history = useHistory();
  return (
    <Grid container justify="center" direction="row" alignItems="flex-start">
      <Grid item xs={12} sm={12} xl={12}>
        <Button
          startIcon={<ArrowBackIosIcon />}
          onClick={() => {
            history.push("/search");
          }}
        >
          Back
        </Button>
      </Grid>
      <Grid item xs={10} sm={12} md={8} xl={6}>
        <Molecule />
      </Grid>
    </Grid>
  );
}

function App() {
  return (
    <div className="App">
      <ThemeProvider theme={theme}>
        <Router>
          <TopBar />
          <Grid container direction="row" justify="center" alignItems="stretch">
            <Grid item xs={12} sm={12} md={10} xl={10}>
              <Suspense fallback={<CircularProgress />}>
                <Switch>
                  <Route exact path="/">
                    <Dashboard />
                  </Route>
                  <Route path="/search">
                    <Search />
                  </Route>
                  <Route path="/about">
                    <About />
                  </Route>
                  <Route path="/download">
                    <Download />
                  </Route>
                  <Route path="/contact">
                    <Contact />
                  </Route>
                  <Route path="/molecules">
                    <h1>Molecules</h1>
                  </Route>
                  <Route path="/dashboard">
                    <Dashboard />
                  </Route>
                  <Route path="/conformer/:conformerId">
                    <Conformer />
                  </Route>
                  <Route path="/molecule/:moleculeId">
                    <MoleculePage />
                  </Route>
                </Switch>
              </Suspense>
            </Grid>
          </Grid>
        </Router>
      </ThemeProvider>
    </div>
  );
}

export default App;
