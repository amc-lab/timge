import { createSlice } from "@reduxjs/toolkit";

const initialState = {
  isLoading: false,
};

const siteSlice = createSlice({
  name: "site",
  initialState,
  reducers: {
    setLoading: (state, action) => {
      state.isLoading = action.payload;
    }
  }
});

export const { setLoading } = siteSlice.actions;
export const selectIsLoading = (state) => state.site.isLoading;
export default siteSlice.reducer;