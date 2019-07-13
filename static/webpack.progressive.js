const merge = require('webpack-merge');
const prod = require('./webpack.prod.js');
const WorkboxPlugin = require('workbox-webpack-plugin');

module.exports = merge(prod, {
    plugins: [
        new WorkboxPlugin.GenerateSW({
            clientsClaim: true,
            skipWaiting: true
        })
    ]
});