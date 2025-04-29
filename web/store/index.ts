import { configureStore } from '@reduxjs/toolkit';
import spaceReducer from './features/space/spaceSlice';
import fileReducer from './features/files/fileSlice';
import siteReducer from './features/site/siteSlice';

export const store = configureStore({
  reducer: {
    space: spaceReducer,
    files: fileReducer,
    site: siteReducer,
  },
});

export type RootState = ReturnType<typeof store.getState>;
export type AppDispatch = typeof store.dispatch;
