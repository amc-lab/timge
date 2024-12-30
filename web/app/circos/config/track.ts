export enum TrackType {
    Karyotype = "Karyotype",
    Chord = "Chord",
}

export interface Track {
    trackType: TrackType;
    data: any;
    config: any;
}