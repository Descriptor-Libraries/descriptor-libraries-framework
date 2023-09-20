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
    forDFT
) {
    return { name, description, forXTB, forDFT }
};

const rows = [
    createData("boltz", "Boltzmann-weighted average of all conformers' properties (T=298.15 K", true, true),
    createData("max", "highest value of a property of any conformer", true, true),
    createData("min", "lowest value of a property of any conformer", true, true),
    createData("std", "standard deviation of the value across all conformers", true, false),
    createData("vburminconf", "property value of the conformer with the smallest buried volume", false, true),
    createData("delta", "difference between the maximum and minimum property values", false, true),
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
                            <TableCell>For XTB</TableCell>
                            <TableCell>For DFT</TableCell>
                        </TableRow>
                    </TableHead>
                    <TableBody>
                        {rows.map((row) => (
                            <TableRow key={row.name}>
                                <TableCell>{row.name}</TableCell>
                                <TableCell>{row.description}</TableCell>
                                <TableCell>{row.forXTB ? '✔️' : ''}</TableCell>
                                <TableCell>{row.forDFT ? '✔️' : ''}</TableCell>
                            </TableRow>
                        ))}
                    </TableBody>
                </Table>
            </TableContainer>
        </Container>
    );
};

export default DataTable;

