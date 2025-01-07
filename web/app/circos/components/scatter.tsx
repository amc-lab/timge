import { useEffect } from "react";
import * as d3 from "d3";
import { ScatterData, ScatterConfig, GlobalConfig } from "../../types/genomes";

interface ScatterProps {
    data: {
        points: Array<ScatterData>;
        globalConfig: GlobalConfig;
        divRef: any;
    };
    config: ScatterConfig;
    segments: Array<any>;
    idx: number;
}

const Scatter = ({ data, config, segments, idx }: ScatterProps) => {
    const canvasRef = data.divRef;

    useEffect(() => {
        let svg = d3.select(canvasRef.current).select("svg");

        if (svg.empty()) {
            svg = d3.select(canvasRef.current)
                .append("svg")
                .attr("width", data.globalConfig.canvasWidth)
                .attr("height", data.globalConfig.canvasHeight);
        }

        const uniqueGroupClass = `group-${idx}`;

        svg.selectAll(`g.${uniqueGroupClass}`).remove();

        const group = svg.append("g")
            .attr("class", uniqueGroupClass)
            .attr("transform", `translate(${data.globalConfig.canvasWidth / 2}, ${data.globalConfig.canvasHeight / 2})`);

        const innerRadius = config.innerRadius;
        const outerRadius = config.innerRadius + 0.8 * config.trackWidth;

        const pointRadiusScale = d3.scaleLinear()
            .domain([0, d3.max(data.points, (d) => d.value) || 100])
            .range([2, 6]);

        group.selectAll("circle.scatter")
            .data(data.points)
            .enter()
            .append("circle")
            .attr("class", "scatter")
            .attr("cx", (d) => {
                const segment = segments.find(seg => seg.chromosome === d.chromosome);
                if (!segment) return null;

                const angle = segment.startAngle + (d.position / segment.length) * (segment.endAngle - segment.startAngle);
                return (innerRadius + (outerRadius - innerRadius) / 2) * Math.cos(angle - Math.PI / 2);
            })
            .attr("cy", (d) => {
                const segment = segments.find(seg => seg.chromosome === d.chromosome);
                if (!segment) return null;

                const angle = segment.startAngle + (d.position / segment.length) * (segment.endAngle - segment.startAngle);
                return (innerRadius + (outerRadius - innerRadius) / 2) * Math.sin(angle - Math.PI / 2);
            })
            .attr("r", (d) => pointRadiusScale(d.value / 100))
            .attr("fill", (d) => {
                const segment = segments.find(seg => seg.chromosome === d.chromosome);
                return segment ? segment.colour : "black";
            })
            .attr("stroke", "black");

    }, [canvasRef, config, data, segments]);

    return null;
};

export default Scatter;
