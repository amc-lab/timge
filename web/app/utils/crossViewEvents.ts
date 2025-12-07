export type CrossViewEventMap = {
  GENERATE_HEATMAP: {
    reference: string;
    track: string;
    segmentA: any;
    segmentB: any;
    resolution: number;
  };
  DELETE_VIEW: { viewId: string };
  ADD_CONNECTION: { source: string; target: string };
  REMOVE_CONNECTION: { viewId: string };
  PROPAGATE_DEPENDENCIES: { viewId: string; dependencies: any };
  PROPAGATE_LOCI: { viewId: string; loci: any };
  UPDATE_VIEW: { index: number; updated: any };
  SET_DIFF_FORM_OPEN: { open: boolean };
  GENERATE_RNAFOLD: {
    reference: string;
    segmentA: string;
    segmentB: string;
    segmentAStart: number;
    segmentAEnd: number;
    segmentBStart: number;
    segmentBEnd: number;
  };
};

type EventKey = keyof CrossViewEventMap;
type Handler<K extends EventKey> = (payload: CrossViewEventMap[K]) => void;

const listeners: { [K in EventKey]?: Set<Handler<K>> } = {};

export const subscribeCrossViewEvent = <K extends EventKey>(
  type: K,
  handler: Handler<K>,
) => {
  if (!listeners[type]) {
    listeners[type] = new Set();
  }
  listeners[type]!.add(handler);
  return () => {
    listeners[type]?.delete(handler);
  };
};

export const publishCrossViewEvent = <K extends EventKey>(
  type: K,
  payload: CrossViewEventMap[K],
) => {
  listeners[type]?.forEach((handler) => {
    handler(payload);
  });
};
