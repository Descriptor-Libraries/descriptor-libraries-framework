import React from "react";

function CDKDepict(props) {
  const [svg, setSVG] = React.useState("");

  const fetchData = () => {
    let url = `/depict/cow/svg?smi=${encodeURIComponent(props.smiles)}`;
    url += `&w=${props.width}`;
    url += `&h=${props.height}`;
    url += `&zoom=${props.zoom}`;
    url += `&abbr=${props.abbr}`;
    url += `&hdisp=${props.hdisp}`;
    url += `&showtitle=${props.showtitle}`;
    url += `&abbr=${props.annotate}`;
    fetch(url)
      .then((r) => r.text())
      .then((r) => {
        r = r.replace(/mm/g, "");
        setSVG(r);
      });
  };

  React.useEffect(fetchData, [props.smiles]);

  return <div dangerouslySetInnerHTML={{ __html: svg }} />;
}

CDKDepict.defaultProps = {
  smiles: "COC1=C(C(OC)=CC=C1)C2=C(P(C3CCCCC3)C4CCCCC4)C=CC=C2",
  width: 80,
  height: 80,
  zoom: 1.6,
  abbr: "off",
  showtitle: false,
  hdisp: "bridgehead",
  annotate: "none",
};

export default CDKDepict;
