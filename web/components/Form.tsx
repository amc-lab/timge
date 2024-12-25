import React from "react";
import { AssemblyConfig, ChordConfig, GlobalConfig } from "@/app/types/genomes";
import Card from '@mui/joy/Card';
import Tabs from '@mui/joy/Tabs';
import Tab from '@mui/joy/Tab';
import TabList from '@mui/joy/TabList';
import TabPanel from '@mui/joy/TabPanel';
import { AssemblyForm } from "./AssemblyForm";
import { ChordForm } from "./ChordForm";
import { GlobalForm } from "./GlobalForm";

interface FormProps {
    handleAssemblyConfigUpdate: (newConfig: AssemblyConfig) => void;
    defaultAssemblyConfig: AssemblyConfig;
    handleChordConfigUpdate: (newConfig: ChordConfig) => void;
    defaultChordConfig: ChordConfig;
    handleGlobalConfigUpdate: (newConfig: GlobalConfig) => void;
    defaultGlobalConfig: GlobalConfig;
}

const Form: React.FC<FormProps> = ({ handleAssemblyConfigUpdate, defaultAssemblyConfig, handleChordConfigUpdate, defaultChordConfig, handleGlobalConfigUpdate, defaultGlobalConfig }) => {
    return (
        <div>
            <Card
                sx={{
                    minWidth: 500,
                    width: '100%',
                }}
            >
                <Tabs
                    orientation="horizontal"
                    size="md"
                >
                    <TabList>
                        <Tab
                            variant="plain"
                            color="neutral"
                            disableIndicator>
                            Assembly
                        </Tab>
                        <Tab
                            variant="plain"
                            color="neutral"
                            disableIndicator>
                            Chords
                        </Tab>
                        <Tab
                            variant="plain"
                            color="neutral"
                            disableIndicator>
                            Global
                        </Tab>
                    </TabList>
                    <TabPanel value={0}>
                        <AssemblyForm onUpdate={handleAssemblyConfigUpdate}
                            defaultConfig={defaultAssemblyConfig}
                        />
                    </TabPanel>
                    <TabPanel value={1}>
                        <ChordForm onUpdate={handleChordConfigUpdate}
                            defaultConfig={defaultChordConfig}
                        />
                    </TabPanel>
                    <TabPanel value={2}>
                        <GlobalForm onUpdate={handleGlobalConfigUpdate}
                            defaultConfig={defaultGlobalConfig}
                        />
                    </TabPanel>
                </Tabs>
            </Card>
        </div>
    )
}

export default Form;