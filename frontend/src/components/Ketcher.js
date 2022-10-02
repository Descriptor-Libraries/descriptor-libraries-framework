import React from 'react'

import { StandaloneStructServiceProvider } from 'ketcher-standalone'
import { Editor } from 'ketcher-react'

import "ketcher-react/dist/index.css";


const structServiceProvider = new StandaloneStructServiceProvider()

class KetcherSketcher extends React.Component {
    render () {
        return (
            <div>
                <Editor
                staticResourcesUrl={process.env.PUBLIC_URL}
                structServiceProvider={structServiceProvider}
                />
            </div>
    )
    }
}

export default KetcherSketcher;