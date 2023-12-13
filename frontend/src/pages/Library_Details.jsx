import {
    Container, 
    Typography, 
    List, 
    ListItem, 
    Link,
    Box,
} from '@mui/material';

import DataTable from '../components/DFTTable';
import ResponsiveYouTube from '../components/YouTube';

import KrakenWorkflow from '../images/kraken-workflow.png';

function Library_Details() {

    return  (
    <Container maxWidth="xl" sx={{marginBottom: "30px"}}>

        <Typography variant="h2" component="h1" textAlign="center">Library Details</Typography>

        <Box sx={{p:4, textAlign:'center'}}>
            <img src={KrakenWorkflow} alt="Kraken Workflow" style={{width:"80%"}}/>
        </Box>
        
        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique",  marginBottom:"15px"}} >Conformational Searching</Typography>
        <Typography variant="body1" paragraph>
        Initial ligand geometries were generated using SMILES strings and converted to free ligands and [LNi(CO)<sub>3</sub>] complexes using RDKit, OpenBabel or Molconvert. 
        These guess geometries were optimized at the GFN2-xTB level of theory (xTB v6.2.2). 
        Optimized geometries were subjected to a conformational search using CREST (v2.8) at the GFN2-xTB level with toluene solvation (GBSA model). 
        For ligands containing ferrocene, conformational searches were performed at the GFN-FF level to avoid structural changes.
        </Typography>
        
        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique",  marginBottom:"15px"}} >xTB Descriptors</Typography>
        <Typography variant="body1" paragraph>
            Molecular descriptors at the xTB level for the full conformational ensembles from CREST of the free ligands and the [NiCO<sub>3</sub>]-bound complexes were collected using MORFEUS.
        </Typography>

        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique", marginTop:"30px", marginBottom:"15px"}}>Selection of conformers for DFT computations</Typography>
        <Typography variant="body1" paragraph>
            Conformers from both sets (free ligand and Ni complex) were selected based on the following two criteria:

            <List sx = {{listStyleType: 'disc', pl: 2,'& .MuiListItem-root': {display: 'list-item',},}}>
            
                <ListItem>
                    Conformers that minimize or maximize any of the following xTB-level steric properties in any of the two conformer sets: 
                    B1, B5, lval, far_vbur, far_vtot, max_delta_qvbur, max_delta_qvtot, near_vbur, near_vtot, ovbur_max, ovbur_min, ovtot_max, 
                    ovtot_min, pyr_val, qvbur_max, qvbur_min, qvtot_max, qvtot_min, vbur.
                </ListItem>
                
                <ListItem>
                Up to 20 conformers within 3 kcal/mol relative energy in the <i>free ligand conformer set.</i> If more than 20 conformers are in that range, the selection was made by RMSD pruning (using PyDP4). 
                This enables structurally diverse selection of conformers in the relevant energy window.
                </ListItem>
            </List>
        </Typography>

        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique", marginTop:"30px", marginBottom:"15px"}}>DFT computations</Typography>

        <Typography variant="body1" paragraph>
                Prior to DFT computations, the [Ni(CO)<sub>3</sub>]-fragment was removed from the Ni complex conformer set to obtain free-ligand initial geometries. 
                All DFT optimizations (Gaussian 16, rev C.01) were performed at the PBE-D3(BJ)/6-31+G(d,p) level of theory. 
                The corresponding geometries were used for a series of single-point energy calculations at the PBE0-D3(BJ)/def2-TZVP and PBE0-D3(BJ)/def2-TZVP/SMD(CHCl<sub>3</sub>) levels of theory. 
                Additional single-points were also run for the radical cations and radical anions from the optimized geometry of the neutral free ligand.
        </Typography>

        <Typography variant="body1" paragraph>
            From the DFT calculations, steric, electronic, or full molecule/interaction-type descriptors were collected for each conformer. 
            The range of properties across the conformers of a single ligand was treated by using up to five condensed measures for each of the properties.
        </Typography>


        <DataTable />
        
        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique",  marginBottom:"15px"}} >ML Descriptor Prediction</Typography>
        
        <Typography variant="body1" paragraph>
            The machine learning portion of the workflow aims to expand the total number of monophosphine ligands in the library from ~1500 (which is a fraction of the total possible space) to {">"} 300,000 ligands. 
            The computational workflow used to calculate the set of ~1500 is too intensive to calculate all possible monophosphines, so descriptors were predicted using machine learning methods.
        </Typography>

        <Typography variant="body1" paragraph>
        New descriptors were predicted using the sum of contributions of each phosphorus substituent using the “Bag of Substituents” approach. 
        These contributions are made up of 576 unique substituents from the ~1500 set of ligands.
        </Typography>

        <Typography variant="body1" paragraph>
            For more details on the machine learning models used to predict descriptors, see the <Link target="_blank" href="https://pubs.acs.org/doi/10.1021/jacs.1c09718">original publication supporting information.</Link>
        </Typography>

        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique",  marginBottom:"15px"}} >Video Tutorial</Typography>
        A video tutorial of Kraken which explains relevant background information and the computational workflow, and provides a walkthrough of this website:
        <Box sx={{p:4, textAlign:'center'}}>
            <ResponsiveYouTube videoId="ApWO7OSvUTk?si=WldA8dwT8Y4WS0_N" />
        </Box>

        <Typography variant="body1" paragraph>
            An extended explanation of the computational workflow used to build the monophosphine library,
            as well as details on the collected and predicted descriptors can be found in <Link target="_blank" href="https://pubs.acs.org/doi/10.1021/jacs.1c09718">original publication supporting information.</Link>
        </Typography>

        <Typography variant="h5" component="h2"  sx={{fontStyle:"oblique",  marginBottom:"15px"}} >Citation</Typography>
        <Typography variant="body1" paragraph>
            Gensch, T.; dos Passos Gomes, G.; Friederich, P.; Peters, E.; Gaudin, T.; Pollice, R.; Jorner, K.; Nigam, A.; Lindner-D’Addario, M.; Sigman, M. S.; Aspuru-Guzik, A.  A Comprehensive Discovery Platform for Organophosphorus Ligands for Catalysis. <span style={{ fontStyle: 'italic' }}>J. Am. Chem. Soc.</span> <span style={{ fontWeight: 'bold' }}>2022</span>, <span style={{ fontStyle: 'italic' }}>144</span>, 3, 1205–1217. DOI: 10.1021/jacs.1c09718
        </Typography>

    </Container>
    )

}

export default Library_Details;