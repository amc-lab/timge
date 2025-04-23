import { getDefaultState } from "../../store/features/space/spaceSlice";

export const STATE_KEY = "timge-state";

export const loadState = () => {
    console.log("Loading state from localStorage");
  
    if (typeof window === "undefined") return getDefaultState();
    console.log("window is defined");
  
    const savedState = localStorage.getItem(STATE_KEY);
  
    if (savedState) {
      try {
        const parsedState = JSON.parse(savedState);
        const defaultState = getDefaultState();
  
        const mergedState = deepMergeDefaults(defaultState, parsedState);
        saveToLocalStorage(mergedState);
        console.log("Merged state:", mergedState);
        return mergedState;
      } catch (error) {
        console.error("Error parsing state from localStorage:", error);
        return getDefaultState();
      }
    }
  
    return getDefaultState();
  };
  
  export const deepMergeDefaults = (defaults: any, target: any): any => {
    if (typeof defaults !== "object" || defaults === null) return target;
    if (typeof target !== "object" || target === null) return defaults;
  
    const result: any = Array.isArray(defaults) ? [...target] : { ...target };
  
    for (const key in defaults) {
      if (!(key in target)) {
        result[key] = defaults[key];
      } else if (
        typeof defaults[key] === "object" &&
        typeof target[key] === "object" &&
        !Array.isArray(defaults[key]) &&
        !Array.isArray(target[key])
      ) {
        result[key] = deepMergeDefaults(defaults[key], target[key]);
      } else {
        result[key] = target[key];
      }
    }
  
    return result;
  };

export const saveToLocalStorage = (newState: any) => {
    if (typeof window !== "undefined") {
        localStorage.setItem(STATE_KEY, JSON.stringify(newState));
    }
}

export const exportState = () => {
    const state = loadState();
    const blob = new Blob([JSON.stringify(state)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "timge-state.json";
    a.click();
    URL.revokeObjectURL(url);
}

