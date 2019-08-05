import React from 'react'
import App from "./App";

const ImportFromFileBodyComponent = () => {
    let fileReader;

    const handleFileRead = (e) => {
        const content = fileReader.result;
        console.log(content);
    };

    const handleFileChosen = (file) => {
        fileReader = new FileReader();
        fileReader.onloadend = handleFileRead;
        fileReader.readAsText(file);
    };

    return <div className={"upload-file"}>
        <input type={'file'} id={'file'} className={'input-file'} accept={".csv"} onChange={
            e => handleFileChosen(e.target.files[0])
        }
               />
    </div>
};

export default ImportFromFileBodyComponent;