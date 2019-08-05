import React, { Component } from 'react'

class App extends Component { // eslint-disable-line react/prefer-stateless-function
    render() {
        return (
                <section>
                        {this.props.children}
                </section>
        );
    }
}

export default App;