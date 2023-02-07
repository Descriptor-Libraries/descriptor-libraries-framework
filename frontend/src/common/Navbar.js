import React from 'react';
import AppBar from '@mui/material/AppBar';
import Box from '@mui/material/Box';
import Toolbar from '@mui/material/Toolbar';
import IconButton from '@mui/material/IconButton';
import Typography from '@mui/material/Typography';
import Menu from '@mui/material/Menu';
import MenuIcon from '@mui/icons-material/Menu';
import Container from '@mui/material/Container';
import Button from '@mui/material/Button';
import MenuItem from '@mui/material/MenuItem';

import Home from '../pages/Home'
import About from '../pages/About'
import Contact from '../pages/Contact'
import Download from '../pages/Download'
import Search  from '../pages/SearchHook';

import {
  BrowserRouter as Router, 
  Routes,
  Route,
  Link
} from 'react-router-dom'

import OriginalKraken from './OriginalKraken';

const pages ={
              'Home': <Home />, 
              'About': <About />, 
              'Search': <Search />, 
              'Download': <Download />,
              'Contact': <Contact />
            };

const ResponsiveAppBar = () => {
  const [anchorElNav, setAnchorElNav] = React.useState(null);

  const handleOpenNavMenu = (event) => {
    setAnchorElNav(event.currentTarget);
  };

  const handleCloseNavMenu = () => {
    setAnchorElNav(null);
  };

  return (
    <Router>
      <AppBar position="static" style={{backgroundColor: "#ed1c24"}} sx={{mb: 2}} >
        <Container maxWidth="xl">
          <Toolbar disableGutters>
            <OriginalKraken sx={{ display: { xs: 'none', md: 'flex' }, mr: 1, fontSize: '60px' }} />
            <Typography
              variant="h6"
              noWrap
              component="a"
              href="/"
              sx={{
                mr: 2,
                display: { xs: 'none', md: 'flex' },
                fontFamily: 'monospace',
                fontWeight: 700,
                letterSpacing: '.3rem',
                color: 'inherit',
                textDecoration: 'none',
              }}
            >
              kraken
            </Typography>

            <Box sx={{ flexGrow: 1, display: { xs: 'flex', md: 'none' } }}>
              <IconButton
                size="large"
                aria-label="account of current user"
                aria-controls="menu-appbar"
                aria-haspopup="true"
                onClick={handleOpenNavMenu}
                color="inherit"
              >
                <MenuIcon />
              </IconButton>
              <Menu
                id="menu-appbar"
                anchorEl={anchorElNav}
                anchorOrigin={{
                  vertical: 'bottom',
                  horizontal: 'left',
                }}
                keepMounted
                transformOrigin={{
                  vertical: 'top',
                  horizontal: 'left',
                }}
                open={Boolean(anchorElNav)}
                onClose={handleCloseNavMenu}
                sx={{
                  display: { xs: 'block', md: 'none' },
                }}
              >
                {Object.keys(pages).map((page) => (
                  <MenuItem id={page} key={page} onClick={handleCloseNavMenu}>
                    <Link to={page.toLowerCase()} id={page} className="NavbarLink"><Typography textAlign="center">{page}</Typography></Link>
                  </MenuItem>
                ))}
              </Menu>
            </Box>
            <OriginalKraken sx={{ display: { xs: 'flex', md: 'none' }, mr: 1, fontSize: '60px' }} />
            <Typography
              variant="h5"
              noWrap
              component="a"
              href=""
              sx={{
                mr: 2,
                display: { xs: 'flex', md: 'none' },
                flexGrow: 1,
                fontFamily: 'monospace',
                fontWeight: 700,
                letterSpacing: '.3rem',
                color: 'inherit',
                textDecoration: 'none',
              }}
            >
              kraken
            </Typography>
            <Box sx={{ flexGrow: 1, display: { xs: 'none', md: 'flex' } }}>
              {Object.keys(pages).map((page) => (
                <Link to={page.toLowerCase()} id={page.toLowerCase()} className="NavbarLink" >
                    <Button
                    key={page}
                    onClick={handleCloseNavMenu}
                    sx={{ my: 2, color: 'white', display: 'block' }}
                  >
                    {page}
                  </Button>
                </Link>
              ))}

              <Link to="/docs/" id="documentation" className="NavbarLink" reloadDocument >
                  <Button
                  key="documentation"
                  sx={{ my: 2, color: 'white', display: 'block' }}>
                  Documentation
                </Button>
              </Link>

            </Box>
          </Toolbar>
        </Container>
      </AppBar>
      
      <Routes>
        <Route path='/' element={[pages['Home']]}></Route>
        {Object.keys(pages).map((page) => (
          <Route path={page} element={pages[page]} >
          </Route>
        )) }
      </Routes>

    </Router>
  );
};
export default ResponsiveAppBar;
