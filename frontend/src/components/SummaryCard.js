import { Typography, Grid, Card, CardContent } from "@mui/material";

export default function StatCard({number, caption, size}) {
    /**
     * Creates a card with a number and caption.
     *  
     *  
     * @param {string} number Number to display in the card.
     * @param {string} caption Caption to display in the card.
     * @param {number} size Size of the card.
     *  
     * @return {StatCard} Card with the number and caption.
     **/

    const cardStyle = {
        height: size,
        display: "flex",
        justifyContent: "center",
        alignItems: "center",
      };


    return (
        <Grid item xs={12} md={6} lg={4}>
                <Card style={ cardStyle }>
                    <CardContent>
                    <Typography variant="h5">{ number }</Typography>
                    <Typography variant="subtitle1"> { caption }</Typography>
                    </CardContent>
                </Card>
        </Grid>
    )
}