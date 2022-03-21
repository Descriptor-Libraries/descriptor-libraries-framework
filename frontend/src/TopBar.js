import React from "react";
import AppBar from "@material-ui/core/AppBar";
import Toolbar from "@material-ui/core/Toolbar";
import Box from "@material-ui/core/Box";
import Typography from "@material-ui/core/Typography";
import IconButton from "@material-ui/core/IconButton";
import GetAppIcon from "@material-ui/icons/GetApp";
import SearchIcon from "@material-ui/icons/Search";
import ScatterPlotIcon from "@material-ui/icons/ScatterPlot";
import InfoIcon from "@material-ui/icons/Info";
import MailOutlineIcon from "@material-ui/icons/MailOutline";
import { makeStyles } from "@material-ui/core/styles";
import useMediaQuery from "@material-ui/core/useMediaQuery";
import { NavLink } from "react-router-dom";

const useStyles = makeStyles((theme) => ({
  title: {
    display: "flex",
    flexGrow: 1,
    color: "#fff",
    marginRight: theme.spacing(2),
    textDecoration: "none",
  },
  logo: {
    marginLeft: theme.spacing(2),
  },
  menuButton: {
    marginRight: theme.spacing(1),
  },
  link: {
    color: "#fff",
  },
  activeLink: {
    color: theme.mtYellow,
  },
}));

export default function TopBar(props) {
  let title = "KRKN";
  if (useMediaQuery("(min-width:400px)")) {
    title = "KRAKEN";
  }
  const theme = useStyles();
  return (
    <AppBar position="static">
      <Toolbar>
        <NavLink to="/" className={theme.title}>
          <img src="/kraken_32x32.png" alt="Kraken logo" />
          <Typography className={theme.title} variant="h6">
            {title}
          </Typography>
        </NavLink>
        <div>
          <NavLink to="/dashboard" className={theme.link}>
            <IconButton className={theme.menuButton} color="inherit" edge="end">
              <ScatterPlotIcon />
              <Box display={{ xs: "none", sm: "none", md: "block" }}>
                <Typography variant="h6">Dashboard</Typography>
              </Box>
            </IconButton>
          </NavLink>
          <NavLink to="/search" className={theme.link}>
            <IconButton className={theme.menuButton} color="inherit" edge="end">
              <SearchIcon />
              <Box display={{ xs: "none", sm: "none", md: "block" }}>
                <Typography variant="h6">Search</Typography>
              </Box>
            </IconButton>
          </NavLink>
          <NavLink to="/about" className={theme.link}>
            <IconButton className={theme.menuButton} color="inherit" edge="end">
              <InfoIcon />
              <Box display={{ xs: "none", sm: "none", md: "block" }}>
                <Typography variant="h6">Info</Typography>
              </Box>
            </IconButton>
          </NavLink>
          <NavLink to="/download" className={theme.link}>
            <IconButton className={theme.menuButton} color="inherit" edge="end">
              <GetAppIcon />
              <Box display={{ xs: "none", sm: "none", md: "block" }}>
                <Typography variant="h6">Download</Typography>
              </Box>
            </IconButton>
          </NavLink>
          <NavLink to="/contact" className={theme.link}>
            <IconButton className={theme.menuButton} color="inherit" edge="end">
              <MailOutlineIcon />
              <Box display={{ xs: "none", sm: "none", md: "block" }}>
                <Typography variant="h6">Contact</Typography>
              </Box>
            </IconButton>
          </NavLink>
        </div>
      </Toolbar>
    </AppBar>
  );
}
