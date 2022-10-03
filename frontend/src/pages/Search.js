import React from 'react';
import { styled } from '@mui/material/styles';
import { TextField, Typography } from "@mui/material";
import Paper from '@mui/material/Paper';
import Grid from '@mui/material/Grid';
import Container from '@mui/material/Container';

import Button from '@mui/material/Button';

const Item = styled(Paper)(({ theme }) => ({
  backgroundColor: theme.palette.mode === 'dark' ? '#1A2027' : '#fff',
  ...theme.typography.body2,
  padding: theme.spacing(1),
  textAlign: 'center',
  color: theme.palette.text.secondary,
}));

class Search extends React.Component {

  // Constructor
  constructor(props) {
    super(props);

    this.state = {
      defaultSearch: 'PC=C',
      searchString: 'PC=C',
      skip: 0,
      limit: 9,
      results: [],
      valid_smiles: true,
    };
  };

  componentDidMount() {
    this.substructureSearch(this.state.defaultSearch)
  }

  dynamicGrid() {
    if (this.state.valid_smiles) {
    return (
    <Container>
      <Grid container spacing={2} sx= {{ mt: 3 }}>
      {
      this.state.results.map((result) => (
        <Grid item xs={12} md={4}>
          <Item>
            <img alt='' src={`data:image/svg+xml;utf8,${encodeURIComponent(result['svg'])}`} />
          </Item> 
        </Grid>
      ))
      }
      
    </Grid>
      <Button variant="contained" style={{backgroundColor: "#ed1c24"}} sx={{ my: 3 }}>Load More</Button>
    </Container>
    )
    }
    else {
      return <Typography sx={{mt:3}}>Invalid SMILES String</Typography>
    }
  }

  substructureSearch(substructure) {
    let encoded = encodeURIComponent(substructure)
    fetch(`/api/v1/molecule/search/?substructure=${encoded}&skip=${this.state.skip}&limit=${this.state.limit}`)
    .then( (response) => {
      if (!response.ok) {
        this.setState({
          valid_smiles: false,
        });
        return [];
      }
      else {
        this.setState({
          valid_smiles: true,
        })
        return response.json()}
      
      })
    .then( (items) => {

      //items.map(item =>  item['svg'] ='');

      this.setState({
        results: items,
        svgs: [],
      })

        items.map( (item) => {
          let encoded = encodeURIComponent(item.smiles);
          let encoded_sub = encodeURIComponent(substructure)
          fetch(`/depict/cow/svg?smi=${encoded}&sma=${encoded_sub}&zoom=1.25&w=50&h=50`)
          .then( response => response.text() )
          .then( (text) => {
            item['svg'] = text;
            this.setState({});
          }
          )
        })


        ;
    } 
   )   
  }

  _handleKeyDown(event) {
    if (event.key === "Enter") {
      this.substructureSearch(this.state.searchString);
    }
  }

  render () {
    return (
    <Container maxWidth="lg">
      <h2>Substructure Search</h2>
      <TextField id="search-outline" 
                label="Enter a SMILES String to Search" 
                variant="outlined"
                 defaultValue= {this.state.defaultSearch} 
                 onChange = { (event) => this.setState({searchString: event.target.value }) }
                 onKeyDown = { (e) => this._handleKeyDown(e) }
                 InputProps={{endAdornment: <Button onClick={ () => { this.substructureSearch(this.state.searchString )} } 
                >
                  Search
                  </Button>}}
                  />
      
      { this.dynamicGrid() }
  </Container>
  
  )
}
}

export default Search;
