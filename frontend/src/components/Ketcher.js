import React, { useEffect, useState } from 'react';
import { Ketcher } from 'ketcher-core';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Editor } from 'ketcher-react';
import Button from '@mui/material/Button';
import Container from '@mui/material/Container';


import "ketcher-react/dist/index.css";


const structServiceProvider = new StandaloneStructServiceProvider()

function KetcherSketcher() {
    const [ketcher, setKetcher] = useState()
    const onClick = () => ketcher.getSmiles().then(result => console.log(result))

        return (
            <Container>
                <Button onClick={onClick}>get value</Button>
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