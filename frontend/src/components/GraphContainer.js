import React, { useState } from 'react';
import Checkbox from '@mui/material/Checkbox';
import FormGroup from '@mui/material/FormGroup';
import FormControlLabel from '@mui/material/FormControlLabel';
import FormControl from '@mui/material/FormControl';
import FormLabel from '@mui/material/FormLabel';
import MyPlot from "../components/Graph"


function GraphContainer() {
  const [checkboxValues, setCheckboxValues] = useState({
    showPCN: true, 
    showPC3: true, 
    showPO3: true, 
    showPN3: true, 
    showPHAL: true 
  });

  const handleCheckboxChange = (event) => {
    setCheckboxValues({
      ...checkboxValues,
      [event.target.name]: !checkboxValues[event.target.name],
    });
  };

  return (
    <div><MyPlot checkboxState={checkboxValues} />
    <FormControl component="fieldset">
      <FormLabel component="legend">Data Visibility</FormLabel>
      <FormGroup aria-label="position" row>
        <FormControlLabel
          name = 'showPC3'
          control={<Checkbox />}
          checked={checkboxValues.showPC3}
          onChange={handleCheckboxChange}
          label="PC3"
          labelPlacement="end"
        />
         <FormControlLabel
          name = 'showPCnRm'
          control={<Checkbox />}
          checked={checkboxValues.showPCnRm}
          onChange={handleCheckboxChange}
          label="P[C]n[R]m"
          labelPlacement="end"
        />
      </FormGroup>
    </FormControl>
    </div>
  );
}

// checkboxState.showPCN !== prevProps.checkboxState.showPCN ||
//       checkboxState.showPC3 !== prevProps.checkboxState.showPC3 ||
//       checkboxState.showPO3 !== prevProps.checkboxState.showPO3 ||
//       checkboxState.showPN3 !== prevProps.checkboxState.showPN3 ||
//       checkboxState.showPHAL !== prevProps.checkboxState.showPHAL

//   return (
//     <div>
//       <Graph checkboxValues={checkboxValues} />
//       <div>
//         <label>
//           <input
//             type="checkbox"
//             name="showPC3"
//             checked={checkboxValues.showPC3}
//             onChange={handleCheckboxChange}
//           />
//           P[C]3
//         </label>
//         <br />
//         <label>
//           <input
//             type="checkbox"
//             name="showPCnRm"
//             checked={checkboxValues.showPCnRm}
//             onChange={handleCheckboxChange}
//           />
//           P[C]3R[m]
//         </label>
//         <br />
//         <label>
//           <input
//             type="checkbox"
//             name="showMarchData"
//             checked={checkboxValues.showPO3}
//             onChange={handleCheckboxChange}
//           />
//           P[O]3
//         </label>
//       </div>
//     </div>
//   );
// }

 export default GraphContainer;
