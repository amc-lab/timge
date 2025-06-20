import { Box, Card, IconButton } from "@mui/joy";
import FolderIcon from "@mui/icons-material/Folder";
import { useAppDispatch } from "@/store/hooks";
import { setDiffStructureFormOpen, setSidebarExpanded } from '@/store/features/space/spaceSlice';

const Sidebar = () => {
    const dispatch = useAppDispatch();
    return (
        <Card
            sx={{
                width: "5%",
                minHeight: "100vh",
                borderRadius: "0",
                padding: "1em",
                overflowY: "auto",
                display: "flex",
                flexDirection: "column",
                alignItems: "center",
                border: "1px solid #ccc", 
                borderTop: "none",
                borderBottom: "none",
            }}
        >
            <Box
                sx={{
                    width: "100%",
                    aspectRatio: "1 / 1",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    marginBottom: "1em",
                }}
            >
                <IconButton
                    variant="plain"
                    color="neutral"
                    sx={{
                        width: "100%",
                        height: "100%",
                        padding: 0,
                    }}
                    onClick={() => dispatch(setSidebarExpanded(true))}
                >
                    <FolderIcon sx={{ width: "60%", height: "60%" }} />
                </IconButton>
            </Box>
            <Box
                sx={{
                    width: "100%",
                    aspectRatio: "1 / 1",
                    display: "flex",
                    alignItems: "center",
                    justifyContent: "center",
                    marginBottom: "1em",
                }}
            >
                <IconButton
                    variant="plain"
                    color="neutral"
                    sx={{
                        width: "100%",
                        height: "100%",
                        padding: 0,
                    }}
                    onClick={() => dispatch(setDiffStructureFormOpen(true))}
                >
                    <FolderIcon sx={{ width: "60%", height: "60%" }} />
                </IconButton>
            </Box>
        </Card>
    );
};

export default Sidebar;
