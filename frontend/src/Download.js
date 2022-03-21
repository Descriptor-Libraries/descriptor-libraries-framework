import Grid from "@material-ui/core/Grid";
import Typography from "@material-ui/core/Typography";

function Download(props) {
  return (
    <div>
      <Grid container direction="row" justify="center" alignItems="stretch">
        <Grid xs={10} sm={10} md={8} lg={6} xl={4} item>
          <Typography variant="h3">Download</Typography>
          <Typography variant="body1" align="justify">
            The data can be browsed through the WebApp. Alternatively, we also
            offer a REST API to query the database, which is documented{" "}
            <a href="http://kraken.cs.toronto.edu/docs">here</a>.
          </Typography>
          <Typography variant="body1" align="justify">
            We will also make the bulk data available shortly.
          </Typography>
        </Grid>
      </Grid>
    </div>
  );
}

export default Download;
