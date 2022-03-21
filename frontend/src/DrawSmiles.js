import React from "react";
import SmilesDrawer from "smiles-drawer";

function DrawSmiles(props) {
  const canvasRef = React.createRef();
  let options = { width: props.width, height: props.height };
  let drawer = new SmilesDrawer.Drawer(options);
  const drawSmiles = (smiles) => {
    if (smiles !== "") {
      SmilesDrawer.parse(smiles, (tree) => {
        drawer.draw(tree, canvasRef.current);
      });
    }
  };

  React.useEffect(() => {
    drawSmiles(props.smiles);
  }, [props.smiles]);

  return (
    <div>
      <canvas
        id="smiles"
        ref={canvasRef}
        width={props.width}
        height={props.height}
      />
    </div>
  );
}

export default DrawSmiles;
