const chordData = [
    {
        "source_chromosome": "os4",
        "source_start": 100000,
        "source_end": 200000,
        "target_chromosome": "os5",
        "target_start": 500000,
        "target_end": 750000
    },
    {
        "source_chromosome": "os11",
        "source_start": 53706,
        "source_end": 119898,
        "target_chromosome": "os12",
        "target_start": 132510,
        "target_end": 318392
    },
    {
        "source_chromosome": "os14",
        "source_start": 890047,
        "source_end": 1176200,
        "target_chromosome": "os9",
        "target_start": 709549,
        "target_end": 839923
    },
    {
        "source_chromosome": "os3",
        "source_start": 533225,
        "source_end": 603174,
        "target_chromosome": "os13",
        "target_start": 589625,
        "target_end": 824058
    },
    {
        "source_chromosome": "os8",
        "source_start": 623161,
        "source_end": 674225,
        "target_chromosome": "os4",
        "target_start": 237627,
        "target_end": 344041
    },
    {
        "source_chromosome": "os1",
        "source_start": 171532,
        "source_end": 213241,
        "target_chromosome": "os4",
        "target_start": 557498,
        "target_end": 632494
    },
    {
        "source_chromosome": "os7",
        "source_start": 611199,
        "source_end": 688226,
        "target_chromosome": "os14",
        "target_start": 1369178,
        "target_end": 1386275
    },
    {
        "source_chromosome": "os5",
        "source_start": 615607,
        "source_end": 711831,
        "target_chromosome": "os3",
        "target_start": 326396,
        "target_end": 406051
    },
    {
        "source_chromosome": "os6",
        "source_start": 337855,
        "source_end": 429067,
        "target_chromosome": "os8",
        "target_start": 241697,
        "target_end": 360185
    },
    {
        "source_chromosome": "os3",
        "source_start": 508749,
        "source_end": 511939,
        "target_chromosome": "os5",
        "target_start": 146236,
        "target_end": 179771
    },
    {
        "source_chromosome": "os4",
        "source_start": 313430,
        "source_end": 414960,
        "target_chromosome": "os14",
        "target_start": 1275583,
        "target_end": 1347285
    },
    {
        "source_chromosome": "os7",
        "source_start": 518785,
        "source_end": 636924,
        "target_chromosome": "os1",
        "target_start": 180067,
        "target_end": 180884
    },
    {
        "source_chromosome": "os4",
        "source_start": 56690,
        "source_end": 83076,
        "target_chromosome": "os14",
        "target_start": 234382,
        "target_end": 325789
    },
    {
        "source_chromosome": "os4",
        "source_start": 515704,
        "source_end": 593275,
        "target_chromosome": "os9",
        "target_start": 274373,
        "target_end": 289908
    },
    {
        "source_chromosome": "os5",
        "source_start": 337479,
        "source_end": 413516,
        "target_chromosome": "os8",
        "target_start": 456526,
        "target_end": 459280
    },
    {
        "source_chromosome": "os9",
        "source_start": 68195,
        "source_end": 140991,
        "target_chromosome": "os4",
        "target_start": 394943,
        "target_end": 412270
    },
    {
        "source_chromosome": "os14",
        "source_start": 184401,
        "source_end": 411744,
        "target_chromosome": "os9",
        "target_start": 483605,
        "target_end": 588244
    },
    {
        "source_chromosome": "os8",
        "source_start": 679191,
        "source_end": 686885,
        "target_chromosome": "os13",
        "target_start": 744636,
        "target_end": 1016986
    },
    {
        "source_chromosome": "os12",
        "source_start": 539417,
        "source_end": 728038,
        "target_chromosome": "os5",
        "target_start": 659441,
        "target_end": 782888
    },
    {
        "source_chromosome": "os5",
        "source_start": 217047,
        "source_end": 227875,
        "target_chromosome": "os8",
        "target_start": 187312,
        "target_end": 297956
    }
]
export const sampleChords = (chordData.map((chord) => ({
        source_chromosome: chord.source_chr,
        source_start: chord.source_start,
        source_end: chord.source_end,
        target_chromosome: chord.target_chr,
        target_start: chord.target_start,
        target_end: chord.target_end,
})));