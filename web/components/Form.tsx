import React from "react";
import { AssemblyConfig, BarConfig, ChordConfig, GlobalConfig } from "@/app/types/genomes";
import Card from '@mui/joy/Card';
import Tabs from '@mui/joy/Tabs';
import Tab from '@mui/joy/Tab';
import TabList from '@mui/joy/TabList';
import TabPanel from '@mui/joy/TabPanel';
import { AssemblyForm } from "./AssemblyForm";
import { ChordForm } from "./ChordForm";
import { GlobalForm } from "./GlobalForm";
import { BarForm } from "./BarForm";
import { Track, TrackType } from "@/app/circos/config/track";
import { defaultGlobalConfig } from "@/app/circos/config/defaultConfigs";

interface FormProps {
    tracks: Array<Track>;
    handleGlobalConfigUpdate: (newConfig: GlobalConfig) => void;
    handleTrackConfigUpdate: (index: number, updatedConfig: any) => void;
}

const Form: React.FC<FormProps> = ({ tracks, handleTrackConfigUpdate, handleGlobalConfigUpdate }) => {
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
                        {
                            tracks.map((track, index) => {
                                return (
                                    <Tab
                                        key={index}
                                        variant="plain"
                                        color="neutral"
                                        disableIndicator>
                                        {track.trackType}
                                    </Tab>
                                )
                            })
                        }
                        <Tab
                            variant="plain"
                            color="neutral"
                            disableIndicator>
                            Global
                        </Tab>
                    </TabList>
                    {
                        tracks.map((track, index) => {
                            if (track.trackType === TrackType.Karyotype) {
                                return (
                                    <TabPanel key={index} value={index}>
                                        <AssemblyForm 
                                            onUpdate={(newConfig: AssemblyConfig) => handleTrackConfigUpdate(index, newConfig)}
                                            defaultConfig={track.config}
                                        />
                                    </TabPanel>
                                )
                            }
                            else if (track.trackType === TrackType.Bar) {
                                return (
                                    <TabPanel key={index} value={index}>
                                        <BarForm 
                                            onUpdate={(newConfig: BarConfig) => handleTrackConfigUpdate(index, newConfig)}
                                            defaultConfig={track.config}
                                        />
                                    </TabPanel>
                                )
                            }
                            else if (track.trackType === TrackType.Chord) {
                                return (
                                    <TabPanel key={index} value={index}>
                                        <ChordForm 
                                            onUpdate={(newConfig: ChordConfig) => handleTrackConfigUpdate(index, newConfig)}
                                            defaultConfig={track.config}
                                        />
                                    </TabPanel>
                                )
                            }
                        })
                    }
                    <TabPanel value={tracks.length}>
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