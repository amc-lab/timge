export class TileCache {
  constructor(limit = 200) {
    this.limit = limit;
    this.map = new Map();
  }

  clear() {
    this.map.clear();
  }

  keys() {
    return Array.from(this.map.keys());
  }

  has(key) {
    return this.map.has(key);
  }

  get(key) {
    const value = this.map.get(key);
    if (!value) return null;

    // Move to end (LRU)
    this.map.delete(key);
    this.map.set(key, value);

    return value;
  }

  set(key, value) {
    if (this.map.has(key)) this.map.delete(key);
    else if (this.map.size >= this.limit) {
      // Evict least recently used
      const firstKey = this.map.keys().next().value;
      this.map.delete(firstKey);
    }
    this.map.set(key, value);
  }
}

export const TILE_SIZE_BINS = 256;

export const buildTileCacheKey = ({
  trackName,
  segmentAId,
  segmentBId,
  resolution,
  tileX,
  tileY,
}) =>
  `${trackName}|${segmentAId}|${segmentBId}|${resolution}|${tileX}|${tileY}`;

export const tileCache = new TileCache(400); // global cache
