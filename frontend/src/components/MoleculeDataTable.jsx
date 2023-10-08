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


async function retrieveData(molecule_id, data_type="ml") {
    try {
        const response = await fetch(`/api/molecules/data/export/${molecule_id}?data_type=${data_type}&return_type=json`);
        const data = await response.json();
        console.log(data_type);
        console.log(data);
        return data;
    } catch (error) {
        console.log(error);
        return null;
    }
}

async function downloadData(molecule_id, data_type) {
    try {
        const response = await fetch(`/api/molecules/data/export/${molecule_id}?data_type=${data_type}&return_type=csv`);
        
        if(response.status === 200) {
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `data_${molecule_id}_${data_type}.csv`;  // you can name the file however you'd like
            a.click();
            window.URL.revokeObjectURL(url);
        } else {
            console.error("Failed to fetch CSV");
        }
    } catch (error) {
        console.log(error);
    }
}


function CustomFooter({ selectedDataType, setSelectedDataType, moleculeID, download }) {
    const handleChange = (event) => {
      setSelectedDataType(event.target.value);
    };
  
    return (
      <GridFooterContainer sx= {{ overflowX: 'auto', whiteSpace: 'nowrap' }}>
        <Select
          value={selectedDataType}
          onChange={handleChange}
          displayEmpty
          sx={{ marginLeft: '8px', marginRight: '16px', height: '40px', verticalAlign: 'middle' }}
        >
          <MenuItem value="ML Data">ML Data</MenuItem>
          <MenuItem value="DFT Data">DFT Data</MenuItem>
          <MenuItem value="XTB Data">XTB Data</MenuItem>
          <MenuItem value="XTB_NI Data">XTB_NI Data</MenuItem>
        </Select>
        <Button disabled={!download} variant="contained" color="primary" sx={{ marginLeft: 'auto', marginRight: '16px', verticalAlign: 'middle' }} onClick={() => { downloadData(moleculeID, dataTypeMapping[selectedDataType]) }} >
            Download CSV
        </Button>
        <GridFooter sx={{
          border: 'none', // To delete double border.
        }} />
      </GridFooterContainer>
    );
  }
  
  function CustomNoRowsOverlay({ selectedDataType}) {
    return (
        <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'center', height: '100%' }}>
            <Typography variant="h6">{selectedDataType} is not available for this molecule.</Typography>
        </div>
    );
}

  

export default function MoleculeDataTable({ molecule_id, initial_data_type }) {
    const [moleculeData, setMoleculeData] = useState(null);
    const [data_type, setDataType] = useState(initial_data_type);
    const [selectedDataType, setSelectedDataType] = useState("ML Data"); // Set the default value to "ML Data"
    const [moleculeID, setMoleculeID] = useState(molecule_id);
    const [download, setDownload] = useState(true); // Allow CSV download

    useEffect(() => {
        async function fetchData() {
            let data = await retrieveData(moleculeID, dataTypeMapping["DFT Data"]);
    
            // If DFT data is empty, default to ML Data
            if (!data || Object.keys(data).length === 0) {
                setSelectedDataType("ML Data");
                data = await retrieveData(moleculeID, dataTypeMapping["ML Data"]);
            } else {
                setSelectedDataType("DFT Data");
            }
            setMoleculeData(data);
        }
    
        fetchData();
    }, [moleculeID]);

    useEffect(() => {
        async function fetchData() {
            const data = await retrieveData(moleculeID, data_type);

            // Check to see if the data is empty, set download to false.
            if (!data) {
                setDownload(false);
            }
            else {
                setDownload(true);
            }
            setMoleculeData(data);
        }

        fetchData();

    }, [data_type]);

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
                components={{Footer: CustomFooter, NoRowsOverlay: CustomNoRowsOverlay}}
                componentsProps={{
                    footer: { selectedDataType, setSelectedDataType, moleculeID, download },
                    noRowsOverlay: { selectedDataType },
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