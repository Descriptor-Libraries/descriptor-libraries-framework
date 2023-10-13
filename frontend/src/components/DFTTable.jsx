import {
    Container, 
    Table,
    TableHead,
    TableContainer,
    TableRow,
    TableBody,
    TableCell,
    Paper,
} from '@mui/material';

function createData(
    name,
    description,
    forXTB,
    forDFT,
    forML,
) {
    return { name, description, forXTB, forDFT, forML }
};

const rows = [
    createData("boltz", "Boltzmann-weighted average of all conformers' properties (T=298.15 K)", true, true, true),
    createData("max", "highest value of a property of any conformer", true, true, true),
    createData("min", "lowest value of a property of any conformer", true, true,true),
    createData("std", "standard deviation of the value across all conformers", true, false, false),
    createData("vburminconf", "property value of the conformer with the smallest buried volume", false, true, true),
    createData("delta", "difference between the maximum and minimum property values", false, true, true),
];

function DataTable() {
    return (
        <Container sx={{marginTop: "20px", marginBottom: "20px"}}>
            <TableContainer component={Paper}>
                <Table>
                    <TableHead>
                        <TableRow>
                            <TableCell>Condensed Properties</TableCell>
                            <TableCell>Description</TableCell>
                            <TableCell>For xTB</TableCell>
                            <TableCell>For DFT</TableCell>
                            <TableCell>For ML</TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {rows.map((row) => (
                            <TableRow key={row.name}>
                                <TableCell>{row.name}</TableCell>
                                <TableCell>{row.description}</TableCell>
                                <TableCell>{row.forXTB ? '✔️' : ''}</TableCell>
                                <TableCell>{row.forDFT ? '✔️' : ''}</TableCell>
                                <TableCell>{row.forML ? '✔️' : ''}</TableCell>
                            </TableRow>
                        ))}
                    </TableBody>
                </Table>
            </TableContainer>
        </Container>
    );
};

export default DataTable;

