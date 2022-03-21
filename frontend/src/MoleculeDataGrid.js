import React from "react";
import Typography from "@material-ui/core/Typography";
import { DataGrid } from "@material-ui/data-grid";
import dft_descriptors from "./dft_descriptors.js";
import xtb_descriptors from "./xtb_descriptors.js";

function MoleculeDataGrid(props) {
  const columns = [
    { field: "method", headerName: "Method", width: 110, type: "string" },
    { field: "name", headerName: "Property", width: 200, type: "string" },
    { field: "descr", headerName: "Description", width: 400, type: "string" },
    { field: "type", headerName: "Type", type: "string" },
    { field: "group", headerName: "Group", hide: true, type: "string" },
    { field: "avg", headerName: "Boltz", width: 90, type: "number" },
    { field: "min", headerName: "Min", width: 90, type: "number" },
    { field: "max", headerName: "Max", width: 90, type: "number" },
    { field: "delta", headerName: "Δ", width: 90, type: "number" },
    { field: "vburmin", headerName: "vburmin", width: 120, type: "number" },
  ];
  const [rows, setRows] = React.useState([]);

  const fmtNum = (n) => {
    if (typeof n == "string") n = +n;
    return n == undefined ? "N/A" : n.toFixed(3);
  };

  React.useEffect(() => {
    let data = [];
    if (props.data.dft_data != undefined) {
      const keys = Object.keys(props.data.dft_data.boltzmann_averaged_data);
      keys.forEach((k, i) =>
        data.push({
          id: i,
          name: k,
          method: "DFT",
          avg: fmtNum(props.data.dft_data.boltzmann_averaged_data[k]),
          min: fmtNum(props.data.dft_data.min_data[k]),
          max: fmtNum(props.data.dft_data.max_data[k]),
          delta: fmtNum(props.data.dft_data.delta_data[k]),
          vburmin: fmtNum(props.data.dft_data.vburminconf_data[k]),
          ...dft_descriptors[k],
        })
      );
    }
    if (props.data.xtb_data != undefined) {
      const size_so_far = data.length;
      const keys = Object.keys(props.data.xtb_data.boltzmann_averaged_data);
      keys.forEach((k, i) =>
        data.push({
          id: i + size_so_far,
          name: k,
          method: "XTB",
          avg: fmtNum(props.data.xtb_data.boltzmann_averaged_data[k]),
          min: fmtNum(props.data.xtb_data.min_data[k]),
          max: fmtNum(props.data.xtb_data.max_data[k]),
          delta: "N/A",
          vburmin: "N/A",
          ...xtb_descriptors[k],
        })
      );
    }
    if (props.data.xtb_ni_data != undefined) {
      const size_so_far = data.length;
      const keys = Object.keys(props.data.xtb_ni_data.boltzmann_averaged_data);
      keys.forEach((k, i) =>
        data.push({
          id: i + size_so_far,
          name: k,
          method: "XTB_NI",
          avg: fmtNum(props.data.xtb_ni_data.boltzmann_averaged_data[k]),
          min: fmtNum(props.data.xtb_ni_data.min_data[k]),
          max: fmtNum(props.data.xtb_ni_data.max_data[k]),
          delta: "N/A",
          vburmin: "N/A",
          ...xtb_descriptors[k],
        })
      );
    }
    if (props.data.ml_data != undefined) {
      const size_so_far = data.length;
      const keys = Object.keys(props.data.ml_data.boltzmann_averaged_data);
      keys.forEach((k, i) =>
        data.push({
          id: i + size_so_far,
          name: k,
          method: "ML",
          avg: fmtNum(props.data.ml_data.boltzmann_averaged_data[k]),
          min: fmtNum(props.data.ml_data.min_data[k]),
          max: fmtNum(props.data.ml_data.max_data[k]),
          delta: "N/A",
          vburmin: "N/A",
          ...dft_descriptors[k],
        })
      );
    }
    setRows(data);
  }, [props.data]);
  return (
    <div style={{ height: 500 }}>
      <Typography variant="h6">Computed data</Typography>
      <DataGrid rows={rows} columns={columns} pageSize={20} />
    </div>
  );
}

export default MoleculeDataGrid;
