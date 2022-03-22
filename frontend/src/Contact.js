import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";

function Contact(props) {
  return (
    <div>
      <Grid container direction="row" justify="center" alignItems="stretch">
        <Grid xs={10} sm={10} md={8} lg={6} xl={4} item>
          <Typography variant="h3">Contact Us</Typography>
          <Typography varuant="body1" align="justify">
            If you have any question or comment, please feel free to reach us at
            the following email: kraken@cs.toronto.edu
          </Typography>
        </Grid>
      </Grid>
    </div>
  );
}

export default Contact;
