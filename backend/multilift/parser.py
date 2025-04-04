from io import StringIO

data = [
    {"chromosome": "os1", "id": 1, "start": 0, "end": 640851, "color": "yellow"},
    {"chromosome": "os2", "id": 2, "start": 0, "end": 947102, "color": "gpos25"},
    {"chromosome": "os3", "id": 3, "start": 0, "end": 1067971, "color": "black"},
    {"chromosome": "os4", "id": 4, "start": 0, "end": 1200490, "color": "blue"},
    {"chromosome": "os5", "id": 5, "start": 0, "end": 1343557, "color": "gpos75"},
    {
        "chromosome": "os6",
        "id": 6,
        "start": 0,
        "end": 1418242,
        "color": "spectral-5-div-4",
    },
    {
        "chromosome": "os7",
        "id": 7,
        "start": 0,
        "end": 1445207,
        "color": "spectral-5-div-4",
    },
    {"chromosome": "os8", "id": 8, "start": 0, "end": 1472805, "color": "gpos100"},
    {"chromosome": "os9", "id": 9, "start": 0, "end": 1541735, "color": "stalk"},
    {"chromosome": "os10", "id": 10, "start": 0, "end": 1687656, "color": "acen"},
    {"chromosome": "os11", "id": 11, "start": 0, "end": 2038340, "color": "green"},
    {"chromosome": "os12", "id": 12, "start": 0, "end": 2271494, "color": "orange"},
    {"chromosome": "os13", "id": 13, "start": 0, "end": 2925236, "color": "blue"},
    {
        "chromosome": "os14",
        "id": 14,
        "start": 0,
        "end": 3291936,
        "color": "spectral-5-div-1",
    },
]


def format_genome(genome, format) -> list[dict]:
    if format == "karyotype" or format == "ring":
        return format_karyotype(genome)
    elif format == "cytoband":
        return format_cytoband(genome)
    elif format == "link":
        return format_link(genome)
    elif format == "bar":
        return format_bar(genome)
    elif format == "heatmap":
        return format_heatmap(genome)
    elif format == "histogram":
        return format_histogram(genome)
    elif format == "scatter":
        return format_scatter(genome)
    elif format == "line":
        return format_line(genome)
    else:
        return []


def format_karyotype(genome) -> list[dict]:
    """
    Args:
    - genome: File object containing the genome sequence

    Genome can be in the following formats:
    - FASTA

    Returns:
    - List of dictionaries containing the karyotype data
    """
    genome_format = genome.name.split(".")[-1]
    if genome_format in ["fasta", "fa"]:
        with StringIO(genome.read().decode("utf-8")) as F:
            lines = F.readlines()
            karyotype = []
            chromosomes = {}
            current_chromosome = None
            current_length = 0

            for line in lines:
                line = line.strip()

                if line.startswith(">"):
                    if current_chromosome:
                        chromosomes[current_chromosome] = current_length
                        karyotype.append(
                            {
                                "chromosome": current_chromosome,
                                "id": len(karyotype) + 1,
                                "start": 0,
                                "end": current_length,
                            }
                        )

                    current_chromosome = line[1:]
                    current_length = 0
                else:
                    current_length += len(line)  # Add to sequence length

            if current_chromosome:
                chromosomes[current_chromosome] = current_length
                karyotype.append(
                    {
                        "chromosome": current_chromosome,
                        "id": len(karyotype) + 1,
                        "start": 0,
                        "end": current_length,
                    }
                )

            return karyotype
    else:
        return []


def format_cytoband(cytoband) -> list[dict]:
    pass


def format_link(link) -> list[dict]:
    """
    Args:
    - link: File object containing the link data

    Returns:
    - List of dictionaries containing the link data

    Link data can be in the following formats:
    - BEDPE
    """

    link_format = link.name.split(".")[-1]
    if link_format in ["bedpe"]:
        with StringIO(link.read().decode("utf-8")) as F:
            lines = F.readlines()
            links = []

            for line in lines:
                line = line.strip()
                data = line.split("\t")
                print(data)
                links.append(
                    {
                        "source_chromosome": data[0],
                        "source_start": data[1],
                        "source_end": data[2],
                        "target_chromosome": data[3],
                        "target_start": data[4],
                        "target_end": data[5],
                        "score": int(data[7]),
                    }
                )

            return links
    else:
        return []


def format_bar(bar) -> list[dict]:
    pass


def format_heatmap(heatmap) -> list[dict]:
    pass


def format_histogram(histogram) -> list[dict]:
    pass


def format_scatter(scatter) -> list[dict]:
    pass


def format_line(line) -> list[dict]:
    line_format = line.name.split(".")[-1]
    if line_format in ["bedgraph"]:
        with StringIO(line.read().decode("utf-8")) as F:
            lines = F.readlines()[1:]
            data = []

            for line in lines:
                line = line.strip()
                data.append(
                    {
                        "chrom": line.split("\t")[0],
                        "chromStart": line.split("\t")[1],
                        "chromEnd": line.split("\t")[2],
                        "value": float(line.split("\t")[3]),
                    }
                )
            return data
    else:
        return []
