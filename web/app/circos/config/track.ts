export enum TrackType {
    Karyotype = "Karyotype",
    Chord = "Chord",
    Bar = "Bar",
    Scatter = "Scatter",
}

export interface Track {
    trackType: TrackType;
    data: any;
    config: any;
}