import { store } from "@/store";
import {
  deleteConnection,
  deleteDependency,
  deleteView,
  setConnection,
  setDependency,
  setDiffStructureFormOpen,
  updateView,
} from "@/store/features/space/spaceSlice";
import { addCustomMapView } from "@/app/utils/viewUtils";
import { subscribeCrossViewEvent } from "./crossViewEvents";
import type { CrossViewEventMap } from "./crossViewEvents";
import { buildApiUrl } from "@/app/config/env";

type Unsubscribe = () => void;

let registered = false;
let unregisterFns: Unsubscribe[] = [];

const handleGenerateHeatmap = (
  payload: CrossViewEventMap["GENERATE_HEATMAP"],
) => {
  const state = store.getState();
  addCustomMapView(store.dispatch, state.space, {
    reference: payload.reference,
    track: payload.track,
    segmentA: payload.segmentA,
    segmentB: payload.segmentB,
    resolution: payload.resolution,
  });
};

const handleDeleteView = (payload: CrossViewEventMap["DELETE_VIEW"]) => {
  store.dispatch(deleteView(payload.viewId));
  store.dispatch(deleteConnection(payload.viewId));
  store.dispatch(deleteDependency(payload.viewId));
};

const handleAddConnection = (
  payload: CrossViewEventMap["ADD_CONNECTION"],
) => {
  const space = store.getState().space;
  const existing = space.connections[payload.source] || [];
  store.dispatch(
    setConnection({
      key: payload.source,
      value: [...existing, payload.target],
    }),
  );
};

const handleRemoveConnection = (
  payload: CrossViewEventMap["REMOVE_CONNECTION"],
) => {
  const space = store.getState().space;
  const dependents = space.connections[payload.viewId];
  if (dependents) {
    dependents.forEach((dependentId) => {
      store.dispatch(setDependency({ key: dependentId, value: [] }));
    });
  }
  store.dispatch(deleteConnection(payload.viewId));
};

const handlePropagateDependencies = (
  payload: CrossViewEventMap["PROPAGATE_DEPENDENCIES"],
) => {
  const space = store.getState().space;
  const dependents = space.connections[payload.viewId];
  if (dependents) {
    dependents.forEach((dependentId) => {
      store.dispatch(
        setDependency({ key: dependentId, value: payload.dependencies }),
      );
    });
  }
};

const handlePropagateLoci = (
  payload: CrossViewEventMap["PROPAGATE_LOCI"],
) => {
  const space = store.getState().space;
  const dependents = space.connections[payload.viewId];
  if (dependents) {
    dependents.forEach((dependentId) => {
      store.dispatch(
        setDependency({
          key: dependentId,
          value: payload.loci,
        }),
      );
    });
  }
};

const handleUpdateView = (payload: CrossViewEventMap["UPDATE_VIEW"]) => {
  store.dispatch(updateView({ index: payload.index, updated: payload.updated }));
};

const handleToggleDiffForm = (
  payload: CrossViewEventMap["SET_DIFF_FORM_OPEN"],
) => {
  store.dispatch(setDiffStructureFormOpen(payload.open));
};

const handleGenerateRNAFold = async (
  payload: CrossViewEventMap["GENERATE_RNAFOLD"],
) => {
  const space = store.getState().space;
  const formData = new URLSearchParams();
  formData.append("uuid", space.uuid);
  formData.append("fasta1", payload.reference);
  formData.append("fasta2", payload.reference);
  formData.append("segment1", payload.segmentA);
  formData.append("segment2", payload.segmentB);
  formData.append(
    "segment1_coords",
    [payload.segmentAStart, payload.segmentAEnd].join(","),
  );
  formData.append(
    "segment2_coords",
    [payload.segmentBStart, payload.segmentBEnd].join(","),
  );

  try {
    const response = await fetch(buildApiUrl("/api/timge/predict_rna_folds/"), {
      method: "POST",
      headers: {
        "Content-Type": "application/x-www-form-urlencoded",
      },
      body: formData.toString(),
    });
    if (!response.ok) {
      throw new Error("Network response was not ok");
    }
    const blob = await response.blob();
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "rnafold_results.gz";
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.URL.revokeObjectURL(url);
  } catch (error) {
    console.error("Error downloading RNAfold results:", error);
  }
};

export const registerCrossViewEventHandlers = () => {
  if (registered) return;
  unregisterFns = [
    subscribeCrossViewEvent("GENERATE_HEATMAP", handleGenerateHeatmap),
    subscribeCrossViewEvent("DELETE_VIEW", handleDeleteView),
    subscribeCrossViewEvent("ADD_CONNECTION", handleAddConnection),
    subscribeCrossViewEvent("REMOVE_CONNECTION", handleRemoveConnection),
    subscribeCrossViewEvent(
      "PROPAGATE_DEPENDENCIES",
      handlePropagateDependencies,
    ),
    subscribeCrossViewEvent("PROPAGATE_LOCI", handlePropagateLoci),
    subscribeCrossViewEvent("UPDATE_VIEW", handleUpdateView),
    subscribeCrossViewEvent("SET_DIFF_FORM_OPEN", handleToggleDiffForm),
    subscribeCrossViewEvent("GENERATE_RNAFOLD", handleGenerateRNAFold),
  ];
  registered = true;
};

export const unregisterCrossViewEventHandlers = () => {
  unregisterFns.forEach((unsub) => unsub());
  unregisterFns = [];
  registered = false;
};
