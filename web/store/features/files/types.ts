export interface File {
    name: string;
    size?: number;
    path: string;
    isDirectory: boolean;
    children?: File[];
}