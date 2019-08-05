let border = {
    "x": 10,
    "y": 20
};

let svgDims = {
    "x": 1000,
    "y": 2000
};

let svgContainer = d3.select("body").append("svg")
    .attr("width", svgDims.x)
    .attr("height", svgDims.y);

let subjectDims = {
    'y': 3,
    'yspacer': 1
};


d3.json("data.json").then(function (data, error) {
    console.log("Data loaded...");
    console.log(data);

    let queryDims = {
        "y": 5,
        "x": data.query.length
    };

    let x = d3.scaleLinear()
        .range([border.x, svgDims.x - border.x])
        .domain([0, queryDims.x]);

    let queryRects = svgContainer.selectAll("query")
        .data([1])
        .enter()
        .append("rect");

    queryRects.attr("class", "query")
        .attr("x", x(0))
        .attr("y", border.y)
        .attr("width", queryDims.x)
        .attr("height", queryDims.y);

    let text = svgContainer.append('text')
    text.text(data.query.name)
        .attr("x", x(queryDims.x/2.0))
        .attr('y', border.y - 5.0)
        .attr("font-family", "sans-serif")
        .attr('text-anchor', 'middle')
        .attr("font-size", "12px");

    let subjectWindow = d3.scaleLinear()
        .range([border.y * 2 + queryDims.y, svgDims.y - border.y])
        .domain([0, svgDims.y]);


    let allsubjects = data.subjects;

    let h = (subjectDims.y + subjectDims.yspacer) * allsubjects.length;

    let y = d3.scaleLinear()
        .range([subjectWindow(0), subjectWindow(h)])
        .domain([0, allsubjects.length]);


    let subjects = [];
    let overOrigin = [];

    for (let i in allsubjects) {
        let s = allsubjects[i];
        if (s.a < s.b) {
            subjects.push([i, s])
        } else {
            overOrigin.push([i, s])
        }
    }

    function handleMouseOver(d, i) {
        console.log("IN!");
        d3.select(this).attr({
            fill: 'orange'
        })
    }

    function handleMouseOut(d, i) {
        console.log("OUT!");
        d3.select(this).attr({
            fill: 'black'
        })
    }

    let subjectRects = svgContainer.selectAll("rects")
        .data(subjects)
        .enter()
        .append('rect')
        .on("mouseover", handleMouseOver)
        .on("mouseout", handleMouseOut);

    subjectRects.attr("class", "subject")
        .attr("height", subjectDims.y)
        .attr('width', function (d) {
            return x(d[1].b) - x(d[1].a)
        })
        .attr("x", function (d) {
            return x(d[1].a)
        })
        .attr("y", function (d) {
            return y(d[0])
        });

    let overOriginRectsLeft = svgContainer.selectAll("rects")
        .data(overOrigin)
        .enter()
        .append('rect');

    overOriginRectsLeft.attr("class", "overOrigin")
        .attr("height", subjectDims.y)
        .attr('width', function (d) {
            return x(d[1].b)
        })
        .attr("x", function (d) {
            return x(0)
        })
        .attr("y", function (d) {
            return y(d[0])
        });

    let overOriginRectsRight = svgContainer.selectAll("rects")
        .data(overOrigin)
        .enter()
        .append('rect');

    overOriginRectsRight.attr("class", "overOrigin")
        .attr("height", subjectDims.y)
        .attr('width', function (d) {
            return x(queryDims.x) - x(d[1].a)
        })
        .attr("x", function (d) {
            return x(d[1].a)
        })
        .attr("y", function (d) {
            return y(d[0])
        });

    let overOriginLines = svgContainer.selectAll("line")
        .data(overOrigin)
        .enter()
        .append('line');

    overOriginLines.attr("class", "overOrigin")
        .attr("x1", function (d) { return x(d[1].b) })
        .attr("x2", function (d) { return x(d[1].a) })
        .attr("y1", function (d) { return y(d[0]) + subjectDims.y * 0.5})
        .attr("y2", function (d) { return y(d[0]) + subjectDims.y * 0.5})
});