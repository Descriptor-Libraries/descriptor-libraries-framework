import * as React from 'react';
import PropTypes from 'prop-types';
import AppBar from '@mui/material/AppBar';
import Box from '@mui/material/Box';
import CssBaseline from '@mui/material/CssBaseline';
import Divider from '@mui/material/Divider';
import Drawer from '@mui/material/Drawer';
import IconButton from '@mui/material/IconButton';
import List from '@mui/material/List';
import ListItem from '@mui/material/ListItem';
import ListItemButton from '@mui/material/ListItemButton';
import ListItemText from '@mui/material/ListItemText';
import MenuIcon from '@mui/icons-material/Menu';
import Toolbar from '@mui/material/Toolbar';
import Typography from '@mui/material/Typography';
import Button from '@mui/material/Button';
import { NavLink, useLocation } from 'react-router-dom';

import { Link } from 'react-router-dom';

import OriginalKraken from './OriginalKraken.js'

const drawerWidth = 240;

function DrawerAppBar(props) {
  const { window, pages } = props;
  const [mobileOpen, setMobileOpen] = React.useState(false);

  const handleDrawerToggle = () => {
    setMobileOpen(prevState => !prevState);
  };

  const location = useLocation();

  const getButtonStyles = (itemPath) => {
    const isHome = (location.pathname === '/' && itemPath === 'home') || location.pathname === `/${itemPath}`;
    return isHome ? { backgroundColor: '#ed1c2477', color: '#fff' } : { color: '#fff' };
  };
  

  const drawer = (
    <Box onClick={handleDrawerToggle} sx={{ textAlign: 'center' }}>
      <Typography variant="h6" sx={{ my: 2 }}>
        kraken
      </Typography>
      <Divider />
      <List>
      {Object.keys(pages).map(item => {
        const itemPath = item.replace(" ", "_").toLowerCase();
        return (
          <ListItem key={item} disablePadding>
            <NavLink 
              to={itemPath} 
              style={{ textDecoration: 'none', color: 'inherit' }}
            >
              <ListItemButton sx={{ textAlign: 'center' }}>
                <ListItemText primary={item} />
              </ListItemButton>
            </NavLink>
          </ListItem>
        );
      })}
      <ListItem key="documentation" disablePadding>
        <Link to="/docs/" id="documentation" className="NavbarLink" reloadDocument style={{ textDecoration: 'none' }} >
          <ListItemButton sx={{ textAlign: 'center' }}>
            <ListItemText primary="Documentation" />
          </ListItemButton>
        </Link>
      </ListItem>
    </List>
  </Box>
  );

  const container = window !== undefined ? () => window().document.body : undefined;

  return (
    <Box sx={{ display: 'flex' }}>
      <CssBaseline />
      <AppBar component="nav">
        <Toolbar>
        <OriginalKraken sx={{ display: { xs: 'none', md: 'flex' }, mr: 1, fontSize: '60px' }} />
          <IconButton
            color="inherit"
            aria-label="open drawer"
            edge="start"
            onClick={handleDrawerToggle}
            sx={{ mr: 2, display: { sm: 'none' } }}
          >
            <MenuIcon />
          </IconButton>
          <Typography
            variant="h6"
            component="div"
            sx={{ flexGrow: 1, display: { xs: 'none', sm: 'block' } }}
          >
            kraken
          </Typography>
          <Box sx={{ display: { xs: 'none', sm: 'block' } }}>
            {Object.keys(pages).map(item => {
              const itemPath = item.replace(" ", "_").toLowerCase();
              return (
                <NavLink
                  to={itemPath}
                  key={item}
                  style={{ textDecoration: 'none' }}
                >
                  <Button sx={getButtonStyles(itemPath)}>
                    <span style={{ textTransform: 'capitalize', fontSize: '16px' }}>{item}</span>
                  </Button>
                </NavLink>
              );
            })}
            <Link to="/docs/" id="documentation" className="NavbarLink" reloadDocument style={{ textDecoration: 'none' }} >
              <Button sx={getButtonStyles('documentation')}>
                <span style={{ textTransform: 'capitalize', fontSize: '16px' }}>Documentation</span>
              </Button>
            </Link>
          </Box>
        </Toolbar>
      </AppBar>
      <nav>
        <Drawer
          container={container}
          variant="temporary"
          open={mobileOpen}
          onClose={handleDrawerToggle}
          ModalProps={{
            keepMounted: true, // Better open performance on mobile.
          }}
          sx={{
            display: { xs: 'block', sm: 'none' },
            '& .MuiDrawer-paper': { boxSizing: 'border-box', width: drawerWidth },
          }}
        >
          {drawer}
        </Drawer>
      </nav>
      <Toolbar />
    </Box>
  );
}

DrawerAppBar.propTypes = {
  window: PropTypes.func, // Injected by the documentation to work in an iframe. You won't need it on your project.
};

export default DrawerAppBar;