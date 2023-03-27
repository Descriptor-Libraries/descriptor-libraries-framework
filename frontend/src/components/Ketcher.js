import React, { useEffect, useState } from 'react';
import { Ketcher } from 'ketcher-core';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Editor } from 'ketcher-react';
import Button from '@mui/material/Button';
import Container from '@mui/material/Container';


import "ketcher-react/dist/index.css";

const structServiceProvider = new StandaloneStructServiceProvider()

function KetcherSketcher() {
    const [ketcher, setKetcher] = useState();
    const [smiles, setSmiles] = useState();
    const [SMARTS, setSMARTS] = useState();

    function onClick(){
        ketcher.getSmiles().then(result => {setSmiles(result);});
    }
    function onClick2(){
        ketcher.getSmarts().then(result => {setSMARTS(result);});
    }
        return (
            <Container>
                <Button onClick={onClick}>get smiles</Button>
                <Button onClick={onClick2}>get smarts</Button>
                <Editor
                    staticResourcesUrl={process.env.PUBLIC_URL}
                    structServiceProvider={structServiceProvider}
                    onInit={(ketcher) => {
                        setKetcher(ketcher)
                    }}
                />
            </Container>
    )       
}

export default KetcherSketcher;