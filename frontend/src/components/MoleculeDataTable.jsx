import React, { useEffect, useState } from 'react';
import { DataGrid, GridFooterContainer, GridFooter } from "@mui/x-data-grid";
import Paper from '@mui/material/Paper';
import Typography from '@mui/material/Typography';
import { Select, MenuItem } from '@mui/material';
import Button from '@mui/material/Button';


const dataTypeMapping = {
    "ML Data": "ml",
    "DFT Data": "dft",
    "XTB Data": "xtb",
    "XTB_NI Data": "xtb_ni"
};


async function retrieveData(molecule_id, data_type) {
    try {
        const response = await fetch(`/api/molecules/data/export/${molecule_id}?data_type=${data_type}&return_type=json`);
        const data = await response.json();
        return data;
    } catch (error) {
        console.log(error);
        return null;
    }
}

function CustomFooter({ selectedDataType, setSelectedDataType }) {
    const handleChange = (event) => {
      setSelectedDataType(event.target.value);
    };
  
    return (
      <GridFooterContainer>
        <Typography sx={{ color: 'gray', display: 'inline-block', verticalAlign: 'middle' }}>
          Data Type:
        </Typography>
        <Select
          value={selectedDataType}
          onChange={handleChange}
          displayEmpty
          sx={{ marginLeft: '8px', marginRight: '16px', display: 'inline-block', verticalAlign: 'middle' }}
        >
          <MenuItem value="ML Data">ML Data</MenuItem>
          <MenuItem value="DFT Data">DFT Data</MenuItem>
          <MenuItem value="XTB Data">XTB Data</MenuItem>
          <MenuItem value="XTB_NI Data">XTB_NI Data</MenuItem>
        </Select>
        <GridFooter sx={{
          border: 'none', // To delete double border.
        }} />
      </GridFooterContainer>
    );
  }
  
  

export default function MoleculeDataTable({ molecule_id, initial_data_type }) {
    const [moleculeData, setMoleculeData] = useState(null);
    const [data_type, setDataType] = useState(initial_data_type);
    const [selectedDataType, setSelectedDataType] = useState("ML Data"); // Set the default value to "ML Data"

    useEffect(() => {
        async function fetchData() {
            let data = await retrieveData(molecule_id, dataTypeMapping["DFT Data"]);
    
            // If DFT data is empty, default to ML Data
            if (!data || Object.keys(data).length === 0) {
                setSelectedDataType("ML Data");
                data = await retrieveData(molecule_id, dataTypeMapping["ML Data"]);
            } else {
                setSelectedDataType("DFT Data");
            }
            setMoleculeData(data);
        }
    
        fetchData();
    }, [molecule_id]);

    useEffect(() => {
        setDataType(dataTypeMapping[selectedDataType]);
        console.log(data_type)
    }, [selectedDataType]);

    const columns = [
        { field: 'property', headerName: 'Property', filterable: true, flex: true },
        { field: 'value', headerName: 'Value', width: 150, filterable: true, headerAlign: 'right', align: 'right', flex: true }
    ];

    const rows = moleculeData ? Object.keys(moleculeData)
        .filter(key => key !== 'smiles' && key !== 'molecule_id')
        .map(key => ({
            id: key,
            property: key,
            value: moleculeData[key],
        })) : [];

    return (
        <Paper elevation={3} style={{ height: 400, width: '100%' }}>
            <DataGrid
                rows={rows}
                columns={columns}
                components={{Footer: CustomFooter}}
                componentsProps={{
                    footer: { selectedDataType, setSelectedDataType }
                }}
                initialState={{
                    pagination: {
                        paginationModel: {
                            pageSize: 25, // This sets the initial page size
                        },
                    },
                }}
                pageSizeOptions={[5, 10, 25, 50, 100]} 
            />
        </Paper>
    );
}