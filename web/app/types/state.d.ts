
import { ViewType } from './viewTypes';

export interface View {
    id: string;
    type: ViewType;
    name: string;
    description: string;
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