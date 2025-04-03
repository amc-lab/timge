"use client"
import Link from 'next/link';
import Image from 'next/image';
import { Box } from '@mui/joy';
import DropdownMenu from './DropdownMenu';
import DropdownButton from './DropdownButton';

interface HeaderButtonProps {
    link: string;
    text: string;
}

const HeaderButton: React.FC<HeaderButtonProps> = ({ link, text }) => {
    return (
        <Box display="flex" alignItems="center" justifyContent="center" height="100%">
            <Link href={link} className="text-white hover:underline">
                {text}
            </Link>
        </Box>
    );
};

interface HeaderProps {
    addLinearGenomeView: () => void;
    addCircosView: () => void;
    addMapView: () => void;
    importTracks: () => void;
    importState: () => void;
    exportState: () => void;
    resetState: () => void;
}

const Header: React.FC<HeaderProps> = ({addLinearGenomeView, addCircosView, addMapView, importTracks, importState, exportState, resetState}) => {
    return (
        <header className="bg-black text-white p-2 flex justify-between items-center">
            <nav>
                <ul className="flex space-x-4 items-center">
                    <li>
                        <DropdownButton label="TIMGE" link="/" />
                    </li>
                    <li>
                        <DropdownMenu 
                            label="File" 
                            items={[
                                { text: "Upload tracks", action: importTracks},
                                { text: "Import state", action: importState },
                                { text: "Export state", action: exportState },
                                { text: "Reset", action: resetState },
                            ]}
                        />
                    </li>
                    <li>
                        <DropdownMenu
                            label="Plot"
                            items={[
                                { text: "Linear View", action: addLinearGenomeView },
                                { text: "Circos View", action: addCircosView },
                                { text: "Map View", action: addMapView },
                            ]}
                        />
                    </li>
                    <li>
                        <DropdownButton label="Multilift" link="/multilift" />
                    </li>
                </ul>
            </nav>
            <nav>
                <ul className="flex space-x-4 items-center">
                    <li>
                        <DropdownButton label="Register" link="/register" />
                    </li>
                    <li>
                        <DropdownButton label="Login" link="/login" />
                    </li>
                </ul>
            </nav>
        </header>
    );
};

export default Header;