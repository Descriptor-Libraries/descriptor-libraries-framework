import React, {useEffect, useState} from 'react';
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TablePagination from '@mui/material/TablePagination';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';
import Paper from '@mui/material/Paper';


async function retrieveData(molecule_id, data_type) {
    try {
        const response = await fetch(`/api/molecules/data/export/${molecule_id}?data_type=${data_type}&return_type=json`);
        const data = await response.json();
        return data;
    } catch (error) {
        debugger;
        console.log(error);
        return null;
    }
}

export default function MoleculeDataTable({molecule_id, initial_data_type}) {
    const [moleculeData, setMoleculeData] = useState(null);
    const [data_type, setDataType] = useState(initial_data_type);
    const [rowsPerPage, setRowsPerPage] = useState(10);
    const [page, setPage] = useState(0);

    useEffect(() => {
        async function fetchData() {
            const data = await retrieveData(molecule_id, data_type);
            setMoleculeData(data);
        }

        fetchData();
    }, [data_type, molecule_id]);
    
  return (
    <TableContainer component={Paper}>
      <Table sx={{ minWidth: 650 }} aria-label="molecule_data_table" stickyHeader>
        <TableHead>
          <TableRow>
            <TableCell>Property</TableCell>
            <TableCell align="right">Value</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
        {moleculeData && 
        Object.keys(moleculeData)
        .filter(key => key !== 'smiles' && key !== 'molecule_id')
        .map((key) => (
            
            <TableRow key={key}>
                <TableCell align="left">{key}</TableCell>
                <TableCell align="right">{moleculeData[key]}</TableCell>
            </TableRow>
            ))};
        </TableBody>
      </Table>
    </TableContainer>
  );
}