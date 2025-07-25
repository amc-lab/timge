export enum TrackType {
  Karyotype = "Karyotype",
  Chord = "Chord",
  Bar = "Bar",
  Ring = "Ring",
  Line = "Line",
  Highlight = "Highlight",
}

export interface Track {
  trackType: TrackType;
  data: any;
  config: any;
  name?: string;
  path?: string;
}
