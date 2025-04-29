import { v4 as uuidv4 } from "uuid";
import { addView, setDependency, setConnection } from "@/store/features/space/spaceSlice";
import { ViewType } from "@/store/features/views/types";

export const addLinearGenomeView = (dispatch, space) => {
    const linearView = {
      uuid: uuidv4(),
      type: ViewType.Linear,
      title: `Linear Genome View ${space.views.length + 1}`,
      description: "Standard linear genome view",
      config: {
        isMinimised: false,
      },
      visible_tracks: [],
    }
    dispatch(addView(linearView));
    dispatch(setDependency({ key: linearView.uuid, value: [] }));
    dispatch(setConnection({ key: linearView.uuid, value: [] }));
}

export const addCircosView = (dispatch, space) => {
    const circosView = {
      uuid: uuidv4(),
      type: ViewType.Circos,
      title: `Circos View ${space.views.length + 1}`,
      description: "Circular genome view",
      visible_tracks: [],
      config: {
        isMinimised: false,
      },
    }
    dispatch(addView(circosView));
    dispatch(setDependency({ key: circosView.uuid, value: [] }));
    dispatch(setConnection({ key: circosView.uuid, value: [] }));
}

export const addMapView = (dispatch, space) => {
    const mapView = {
      uuid: uuidv4(),
      type: ViewType.Map,
      title: `Map View ${space.views.length + 1}`,
      description: "Map view",
      config: {
        reference: "",
        track: "",
        segmentA: "",
        segmentB: "",
        resolution: 5,
        isMinimised: false,
      },
    }
    dispatch(addView(mapView));
    dispatch(setDependency({ key: mapView.uuid, value: [] }));
    dispatch(setConnection({ key: mapView.uuid, value: [] }));
}

export const addCustomMapView = (dispatch, space, mapConfig) => {
    const customMapView = {
      uuid: uuidv4(),
      type: ViewType.Map,
      title: `Map View ${space.views.length + 1}`,
      description: "Map view",
      config: mapConfig,
    }
    dispatch(addView(customMapView));
    dispatch(setDependency({ key: customMapView.uuid, value: [] }));
    dispatch(setConnection({ key: customMapView.uuid, value: [] }));
}