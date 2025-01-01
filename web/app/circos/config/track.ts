export enum TrackType {
    Karyotype = "Karyotype",
    Chord = "Chord",
    Bar = "Bar",
}

export interface Track {
    trackType: TrackType;
    data: any;
    config: any;
}