import { View } from "../views/types";

export interface State {
    title: string;
    uuid: string;
    views: View[];
    dateCreated: string;
    dateModified: string;
    isUserLoggedIn: boolean;
    // mapping of view id to the list of view ids that are connected to it
    connections: Record<string, string[]>;
    // mapping of view id to the list of view ids that are dependent on it
    dependencies: Record<string, string[]>;
    dataFiles: string[];
    config?: Record<string, any>;
}