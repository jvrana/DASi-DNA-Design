import React from 'react'
import {render} from 'react-dom'
import App from './App'
import ImportFromFileBodyComponent from './ImportFile'
import Chart from './TestChart'

render(
    <App>
        <div>
            Hello world!
        </div>
        <ImportFromFileBodyComponent/>
        <Chart></Chart>
    </App>,
    document.getElementById('root')
);