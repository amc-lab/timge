"use client"
import Link from 'next/link';
import Image from 'next/image';
import { Box, Button } from '@mui/joy';
import DropdownMenu from './DropdownMenu';
import DropdownButton from './DropdownButton';
import { useAppDispatch } from "@/store/hooks";
import { setMultiliftFormOpen } from '@/store/features/space/spaceSlice';

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
    const dispatch = useAppDispatch();
    
    return (
<header
    className="bg-black text-white flex justify-between items-center"
    style={{
        height: '3.5rem', // Fixed height
        padding: '0 1rem',
        minWidth: 0,
        flexShrink: 0,
    }}
>
    <nav>
        <ul className="flex items-center space-x-4">
            <li>
                <DropdownButton label="vRISE" link="/" />
            </li>
            <li>
                <DropdownMenu 
                    label="File"
                    items={[
                        { text: "Upload tracks", action: importTracks },
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
                <Button
                    variant="solid"
                    sx={{
                    backgroundColor: "black",
                    color: "white",
                    height: "100%",
                    display: "flex",
                    alignItems: "center",
                    px: 3,
                    "&:hover": { backgroundColor: "#333" },
                    minWidth: "4em"
                    }}
                    onClick={() => dispatch(setMultiliftFormOpen(true))}
                    >
                        Multilift
                </Button>
            </li>
            <li>
                <DropdownButton label="Docs" link="/docs" />
            </li>
        </ul>
    </nav>
    <nav>
        <ul className="flex items-center space-x-4">
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