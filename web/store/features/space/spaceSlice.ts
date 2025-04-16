import { createSlice, PayloadAction } from '@reduxjs/toolkit';
import { State } from './types';
import { v4 as uuidv4 } from 'uuid';
import { loadState } from '@/app/utils/stateUtils';

export const getDefaultState = (): State => ({
  title: 'Untitled',
  uuid: uuidv4(),
  views: [],
  dateCreated: new Date().toISOString(),
  dateModified: new Date().toISOString(),
  isUserLoggedIn: false,
  connections: {},
  dependencies: {},
  dataFiles: [], // Array of file paths or URLs
});

const spaceSlice = createSlice({
  name: 'space',
  initialState: loadState(),
  reducers: {
    setSpace: (state, action: PayloadAction<State>) => action.payload,
    resetSpace: () => getDefaultState(),
    addView: (state, action: PayloadAction<any>) => {
      state.views.push(action.payload);
    },
    deleteView: (state, action: PayloadAction<string>) => {
      state.views = state.views.filter(v => v.uuid !== action.payload);
    },
    updateView: (state, action: PayloadAction<{ index: number, updated: any }>) => {
      state.views[action.payload.index] = {
        ...state.views[action.payload.index],
        ...action.payload.updated,
      };
    },
    addDataFile: (state, action: PayloadAction<string>) => {
      state.dataFiles.push(action.payload);
    },
    deleteDataFile: (state, action: PayloadAction<string>) => {
      state.dataFiles = state.dataFiles.filter(file => file !== action.payload);
    },
    setDataFiles: (state, action: PayloadAction<string[]>) => {
      state.dataFiles = action.payload;
    },
    setDependency: (state, action: PayloadAction<{ key: string, value: any }>) => {
      const { key, value } = action.payload;
      state.dependencies[key] = value;
    },
    deleteDependency: (state, action: PayloadAction<string>) => {
      delete state.dependencies[action.payload];
    },
    setConnection: (state, action: PayloadAction<{ key: string, value: any }>) => {
      const { key, value } = action.payload;
      state.connections[key] = value;
    },
    deleteConnection: (state, action: PayloadAction<string>) => {
      delete state.connections[action.payload];
    },
    setTitle: (state, action: PayloadAction<string>) => {
      state.title = action.payload;
    },
    setViewTitle: (state, action: PayloadAction<{ uuid: string, title: string }>) => {
      const { uuid, title } = action.payload;
      const view = state.views.find(v => v.uuid === uuid);
      if (view) {
        view.title = title;
      }
    },
  },
});

export const {
  setSpace,
  resetSpace,
  addView,
  deleteView,
  updateView,
  addDataFile,
  setDataFiles,
  setDependency,
  setConnection,
  deleteConnection,
  deleteDependency,
  setTitle,
  setViewTitle,
} = spaceSlice.actions;

export const selectSpace = (state: { space: State }) => state.space;

export default spaceSlice.reducer;
