export enum TrackType {
    karotype = "karotype",
    chord = "chord",
}

export interface Track {
    trackType: TrackType;
    data: any;
    config: any;
}