const trimTrailingSlash = (value?: string | null) => {
  if (!value) return undefined;
  const trimmed = value.trim();
  if (!trimmed) return undefined;
  return trimmed.replace(/\/+$/, "");
};

const withDefaultApiHost = trimTrailingSlash(process.env.NEXT_PUBLIC_DJANGO_HOST) ?? "http://127.0.0.1:8000";
const withDefaultFileHost =
  trimTrailingSlash(process.env.NEXT_PUBLIC_FILE_HOST) ?? `${withDefaultApiHost}/uploads`;

export const API_BASE_URL = withDefaultApiHost;
export const FILE_BASE_URL = withDefaultFileHost;

export const buildApiUrl = (path: string) =>
  `${API_BASE_URL}${path.startsWith("/") ? path : `/${path}`}`;

export const buildFileUrl = (path: string) => {
  const normalised = path.replace(/^\/+/, "");
  return `${FILE_BASE_URL}/${normalised}`;
};
