/**
 * Coordinate helpers.
 *
 * The exported data uses EPSG:25830 metres, recentred around the basin centre.
 * In the Three.js scene we work in kilometres:
 *   world X = easting  (east  +)
 *   world Z = -northing (north  -Z, so "up" on screen is north)
 *   world Y = elevation (up)
 */
export const SCALE = 0.001 // scene units per metre

export function sceneX(xMetres) {
  return xMetres * SCALE
}

export function sceneZ(yMetres) {
  return -yMetres * SCALE
}

export function sceneHeight(elevMetres, exaggeration = 1) {
  return elevMetres * exaggeration * SCALE
}

export function unproject(manifest, x, y) {
  return [x + manifest.origin[0], y + manifest.origin[1]]
}

/** Scene-space centre of the recentred data bounds. */
export function sceneCentre(bounds) {
  return [
    (bounds.minX + bounds.maxX) / 2 * SCALE,
    -((bounds.minY + bounds.maxY) / 2) * SCALE,
  ]
}

/** Approximate extent of the data in scene units. */
export function sceneExtent(bounds) {
  return {
    width: (bounds.maxX - bounds.minX) * SCALE,
    depth: (bounds.maxY - bounds.minY) * SCALE,
  }
}
