import { getDefaultState } from "../../store/features/space/spaceSlice";

export const STATE_KEY = "timge-state";

export const loadState = () => {
    console.log("Loading state from localStorage");
    if (typeof window === "undefined") return getDefaultState();
    console.log("window is defined");
    const savedState = localStorage.getItem(STATE_KEY);
    if (savedState) {
        try {
        return JSON.parse(savedState);
        }
        catch (error) {
        console.error("Error parsing state from localStorage:", error);
        return getDefaultState();
        }
    }
    return getDefaultState();
}

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

