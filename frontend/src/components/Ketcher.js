import React, { useEffect, useState } from 'react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import { Editor } from 'ketcher-react';
import Button from '@mui/material/Button';
import Container from '@mui/material/Container';


import "ketcher-react/dist/index.css";


const structServiceProvider = new StandaloneStructServiceProvider()

function KetcherSketcher() {
    const [ketcher, setKetcher] = useState()
    const onClick = () => ketcher.getSmarts().then(result => console.log(result))

        return (
            <Container>
                <Button onClick={onClick}>get value</Button>
                <Editor
                    staticResourcesUrl={process.env.PUBLIC_URL}
                    structServiceProvider={structServiceProvider}
                    onInit={(ketcherInstance) => {
                        setKetcher(ketcherInstance)
                    }}
                />
            </Container>
    )       
}

export default KetcherSketcher;