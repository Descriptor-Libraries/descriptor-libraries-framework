import Grid from "@mui/material/Grid";
import Typography from "@mui/material/Typography";

function About(props) {
  return (
    <div>
      <Grid container direction="row" justify="center" alignItems="stretch">
        <Grid xs={10} sm={10} md={8} lg={6} xl={4} item>
          <Typography variant="h3">About</Typography>
          <Typography variant="subtitle1">
            Kraken stands for ...<b>K</b>olossal vi<b>R</b>tual d<b>A</b>tabase
            for mole<b>K</b>ular d<b>E</b>scriptors of orga<b>N</b>ophosphorus
            ligands.
          </Typography>

          <img src="/kraken.png" alt="Kraken logo" />

          <Typography variant="body1" align="justify">
            With descriptors for 330949 monodentate organophosphorus(III)
            ligands at two levels of theory as well as property estimation
            powered by machine learning, we hope experimentalists,
            theoreticians, and data scientists will find utility in this library
            when designing new ligands for catalysis. This descriptor set
            accounts for conformational flexibility and is continuously
            developed and maintained by the authors of "Mapping the Property
            Space of Monodentate Organophosphorus Ligands for Catalysis",
            preprint on ChemRxiv (doi: 10.26434/chemrxiv.12996665).
          </Typography>

          <Typography variant="body1" align="justify">
            This project is a collaboration between University of Toronto,
            University of Utah, Technische Universität Berlin, Karlsruhe
            Institute of Technology, Vector Institute for Artificial
            Intelligence, Center for Computer Assisted Synthesis, IBM Research,
            and AstraZeneca.
          </Typography>
          
          <a href="https://www.utoronto.ca/">
            <img src="/uoft.png" width="25%" alt-text="University of Toronto" />
          </a>

          <a href="https://www.utah.edu/">
            <img src="/utah.png" width="50%" alt-text="University of Utah" />
          </a>
          <a href="https://www.kit.edu/english/">
            <img
              src="/kit.png"
              width="50%"
              alt-text="Karlsuher Institut für Technologie"
            />
          </a>
          <a href="https://www.tu.berlin/en/">
            <img src="/tu_berlin.png" width="50%" alt-text="TU Berlin" />
          </a>
          <a href="https://vectorinstitute.ai/">
            <img src="/vector.png" width="50%" alt-text="Vector Institute" />
          </a>
          <a href="http://matter.toronto.edu">
            <img src="/matter.png" width="50%" alt-text="Matter Lab" />
          </a>
          <a href="https://chem.utah.edu/directory/sigman/research-group/index.php">
            <img src="/sigman_lab.png" width="30%" alt-text="Sigman Lab" />
          </a>
          <Typography variant="h5">Funding</Typography>
          <Typography variant="body1" align="justify">
            We would like to thank our funders:
            <br />
            G.P.G gratefully acknowledges the Natural Sciences and Engineering
            Research Council of Canada (NSERC) for the Banting Postdoctoral
            Fellowship. P.F. acknowledges funding from the European Union’s
            Horizon 2020 research and innovation programme under the Marie
            Skodowska-Curie grant agreement no. 795206 (MolDesign). K.J. is a
            fellow of the AstraZeneca Postdoc Programme. T.G. thanks the
            Leopoldina Fellowship Programme of the German National Academy of
            Sciences Leopoldina (LPDS 2017−18) and thanks the Center for High
            Performance Computing at the University of Utah for support and
            resources. Funded by the Deutsche Forschungsgemeinschaft (DFG,
            German Research Foundation) under Germany's Excellence Strategy –
            EXC 2008/1 – 390540038. Gefördert durch die Deutsche
            Forschungsgemeinschaft (DFG) im Rahmen der Exzellenzstrategie des
            Bundes und der Länder – EXC 2008/1 – 390540038. We acknowledge the
            Defense Advanced Research Projects Agency (DARPA) under the
            Accelerated Molecular Discovery Program under Cooperative Agreement
            No. HR00111920027 dated August 1, 2019. The content of the
            information presented in this work does not necessarily reflect the
            position or the policy of the Government. A. A.-G. is grateful for
            the support of Anders G. Frøseth, Natural Resources Canada, and the
            Canada 150 Research Chairs program. Additionally, A. A.-G. thanks
            Compute Canada for computational resources. DFT and xtb calculations
            were performed on the niagara supercomputer at the SciNet HPC
            Consortium, which is funded by the Canada Foundation for Innovation;
            the Government of Ontario; Ontario Research Fund - Research
            Excellence; and the University of Toronto. Machine learning models
            were developed on the supercomputer beluga from École de technologie
            supérieure, managed by Calcul Québec and Compute Canada and funded
            by the Canada Foundation for Innovation (CFI), the ministère de
            l’Économie, de la science et de l’innovation du Québec (MESI) and
            the Fonds de recherche du Québec - Nature et technologies (FRQ-NT).
            We are grateful to Dr. Claire Yu for helping with the deployment of
            the web app. M. S. S. thanks the NSF under the CCI Center for
            Computer Assisted Synthesis (CHE-1925607) for support.
          </Typography>
        </Grid>
      </Grid>
    </div>
  );
}

export default About;
