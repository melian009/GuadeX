/*
 * Shared GIS/CSV helpers for the GuadeX data build scripts.
 * ---------------------------------------------------------
 * Used by build-data.mjs and export-crosswalk.mjs so shapefile attribute
 * decoding and feature reading stay in one place.
 */
import shapefile from 'shapefile'

/** Shapefile text is stored as UTF-8 bytes decoded as latin1 -> repair it. */
export function fixText(value) {
  if (typeof value !== 'string') return value
  if (!/[ÃÂ]/.test(value)) return value
  try {
    return Buffer.from(value, 'latin1').toString('utf8')
  } catch {
    return value
  }
}

/** Read a shapefile fully (streaming) as an array of GeoJSON features. */
export async function readFeatures(path) {
  const source = await shapefile.open(path)
  const features = []
  for (;;) {
    const { done, value } = await source.read()
    if (done) break
    if (value) features.push(value)
  }
  return features
}

/** Quote a CSV cell when it contains a delimiter, quote, or newline. */
export function csvCell(value) {
  const s = String(value ?? '')
  return /[",\n\r]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s
}
