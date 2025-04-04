
import { ViewType } from './viewTypes';

export interface View {
    id: string;
    type: ViewType;
    title: string;
    description: string;
    config?: any;
    visible_tracks?: string[];
    uuid?: string;
}

export interface State {
    title: string;
    // dateCreated: string;
    // dateModified: string;
    isUserLoggedIn: boolean;
    UUID: string?;
    // dataFilesRootDir: string;
    dataFiles: string[]; 
    views: View[];
}