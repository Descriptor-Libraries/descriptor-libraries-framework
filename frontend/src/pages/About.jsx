import Typography from '@mui/material/Typography';
import Container from '@mui/material/Container';

function About() {
    return (<>
        <Container maxWidth="lg">
            <h2>About</h2>
            
            <Typography variant="subtitle1" sx={{mb:3}}>
                Kraken stands for ...<b>K</b>olossal vi<b>R</b>tual d<b>A</b>tabase
                for mole<b>K</b>ular d<b>E</b>scriptors of orga<b>N</b>ophosphorus
                ligands.
            </Typography>

            <Typography variant="body1" align="justify" paragraph={true}> 
                With descriptors for 330949 monodentate organophosphorus(III)
                ligands at two levels of theory as well as property estimation
                powered by machine learning, we hope experimentalists,
                theoreticians, and data scientists will find utility in this library
                when designing new ligands for catalysis. This descriptor set
                accounts for conformational flexibility and was created by the authors of "Mapping the Property
                Space of Monodentate Organophosphorus Ligands for Catalysis",
                preprint on ChemRxiv (doi: 10.26434/chemrxiv.12996665).
            </Typography>

            <Typography variant="body1" align="justify" paragraph={true}>
                This project was originally created as a collaboration between University of Toronto,
                University of Utah, Technische Universität Berlin, Karlsruhe
                Institute of Technology, Vector Institute for Artificial
                Intelligence, Center for Computer Assisted Synthesis, IBM Research,
                and AstraZeneca.
            </Typography>

            
            <Typography variant="body1" align="justify" paragraph={true}>
                This project is now maintained by The Molecular Sciences Software Institute.
            </Typography>
        </Container>
    </>) 
}

export default About;
