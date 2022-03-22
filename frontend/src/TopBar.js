import React from "react";
import AppBar from "@mui/material/AppBar";
import Toolbar from "@mui/material/Toolbar";
import Box from "@mui/material/Box";
import Typography from "@mui/material/Typography";
import IconButton from "@mui/material/IconButton";
import GetAppIcon from "@mui/icons-material/GetApp";
import SearchIcon from "@mui/icons-material/Search";
import ScatterPlotIcon from "@mui/icons-material/ScatterPlot";
import InfoIcon from "@mui/icons-material/Info";
import MailOutlineIcon from "@mui/icons-material/MailOutline";
import { makeStyles } from "@mui/material/styles";
import useMediaQuery from "@mui/material/useMediaQuery";
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
