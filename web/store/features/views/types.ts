export enum ViewType {
    Circos = 'circos',
    Linear = 'linear',
    Map = 'map',
}

export interface View {
    // id: string;
    type: ViewType;
    title: string;
    description: string;
    config?: any;
    visible_tracks?: string[];
    uuid: string;
}