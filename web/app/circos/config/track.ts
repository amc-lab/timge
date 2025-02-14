export enum TrackType {
  Karyotype = "Karyotype",
  Chord = "Chord",
  Bar = "Bar",
  Ring = "Ring",
}

export interface Track {
  trackType: TrackType;
  data: any;
  config: any;
}
