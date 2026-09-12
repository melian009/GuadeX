/*
 * GuadeX 3D visualization - data build pipeline
 * ---------------------------------------------
 * Converts the raw Guadalquivir GIS shapefiles and the sampling-site tables
 * into compact, recentred JSON layers consumed by the Three.js viewer.
 *
 * Source CRS: EPSG:25830 (ETRS89 / UTM zone 30N, metres).
 * Output coordinates are recentred around the basin centre and rounded to
 * whole metres, so the browser never sees the large UTM magnitudes.
 *
 * Run with:  npm run data
 */
import { fileURLToPath } from 'node:url'
import { dirname, join } from 'node:path'
import { mkdir, writeFile, stat } from 'node:fs/promises'
import { readFileSync } from 'node:fs'
import shapefile from 'shapefile'

const __dirname = dirname(fileURLToPath(import.meta.url))
const ROOT = join(__dirname, '..')
const DATA = join(ROOT, '..', 'data')
const OUT = join(ROOT, 'public', 'data')

const LAYER_FILES = {
  catchments: join(DATA, 'Version_02-01-2026-Masas', 'Cuencas_masas_agua_4c.shp'),
  gwb: join(DATA, 'Version_02-01-2026-Masas', 'GWB_4C.shp'),
  rivers: join(DATA, 'Version_02-01-2026-Masas', 'SW_Line_4C_.shp'),
  reservoirs: join(DATA, 'Version_02-01-2026-Masas', 'SWB_Lago.shp'),
  transitional: join(DATA, 'Version_02-01-2026-Masas', 'SWB_Transicion.shp'),
  coastal: join(DATA, 'Version_02-01-2026-Masas', 'SWB_Costera.shp'),
  massesPoints: join(DATA, 'GIS', 'capas GIS', 'capa_masas con puntos muestreo', 'masas con puntos muestreo.shp'),
  sites: join(DATA, 'GIS', 'capas GIS', 'capa_puntos muestreo con id masa', 'puntos muestreo con id masa.shp'),
  connectivity: join(DATA, 'ConnectivityUTM.csv'),
}

/* Aggressive simplification tolerances in metres (per layer). */
const TOLERANCE = {
  catchments: 90,
  gwb: 180,
  rivers: 30,
  reservoirs: 35,
  transitional: 40,
  coastal: 60,
  massesPoints: 120,
}

/* ------------------------------------------------------------------ */
/* Helpers                                                            */
/* ------------------------------------------------------------------ */

/** shapefile text is stored as UTF-8 bytes decoded as latin1 -> repair it. */
function fixText(value) {
  if (typeof value !== 'string') return value
  if (!/[ÃÂ]/.test(value)) return value
  try {
    return Buffer.from(value, 'latin1').toString('utf8')
  } catch {
    return value
  }
}

function cleanProps(props) {
  const out = {}
  for (const [k, v] of Object.entries(props)) {
    if (v === null || v === undefined) continue
    if (typeof v === 'string') {
      const t = fixText(v).trim()
      out[k] = t
    } else if (typeof v === 'number') {
      out[k] = Number.isFinite(v) ? Math.round(v * 1e4) / 1e4 : null
    } else {
      out[k] = v
    }
  }
  return out
}

/** Read a shapefile fully as GeoJSON features (streaming, lower memory). */
async function readFeatures(path) {
  const source = await shapefile.open(path)
  const features = []
  for (;;) {
    const { done, value } = await source.read()
    if (done) break
    if (value) features.push(value)
  }
  return features
}

/** Iterative Douglas-Peucker for an open polyline of [x,y] points. */
function simplifyLine(points, tol) {
  if (points.length <= 2) return points.slice()
  const keep = new Uint8Array(points.length)
  keep[0] = 1
  keep[points.length - 1] = 1
  const stack = [[0, points.length - 1]]
  while (stack.length) {
    const [first, last] = stack.pop()
    if (last <= first + 1) continue
    const [x1, y1] = points[first]
    const [x2, y2] = points[last]
    const dx = x2 - x1
    const dy = y2 - y1
    const norm = Math.hypot(dx, dy)
    let maxD = -1
    let idx = -1
    for (let i = first + 1; i < last; i++) {
      const [px, py] = points[i]
      const d = norm === 0
        ? Math.hypot(px - x1, py - y1)
        : Math.abs(dy * px - dx * py + x2 * y1 - y2 * x1) / norm
      if (d > maxD) { maxD = d; idx = i }
    }
    if (maxD > tol && idx !== -1) {
      keep[idx] = 1
      stack.push([first, idx], [idx, last])
    }
  }
  const out = []
  for (let i = 0; i < points.length; i++) if (keep[i]) out.push(points[i])
  return out
}

function ringToPoints(ring) {
  const pts = ring.map((p) => [p[0], p[1]])
  // drop duplicated closing vertex for simplification
  if (pts.length > 1 &&
      pts[0][0] === pts[pts.length - 1][0] &&
      pts[0][1] === pts[pts.length - 1][1]) {
    pts.pop()
  }
  return pts
}

function pointsToRing(pts) {
  if (pts.length < 3) return null
  const ring = pts.slice()
  ring.push([ring[0][0], ring[0][1]])
  return ring
}

function ringArea(pts) {
  let a = 0
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
    a += pts[j][0] * pts[i][1] - pts[i][0] * pts[j][1]
  }
  return Math.abs(a) / 2
}

/* ------------------------------------------------------------------ */
/* Coordinate transform                                               */
/* ------------------------------------------------------------------ */

const state = {
  origin: [0, 0],
  min: [Infinity, Infinity],
  max: [-Infinity, -Infinity],
}

function noteBounds(x, y) {
  if (x < state.min[0]) state.min[0] = x
  if (y < state.min[1]) state.min[1] = y
  if (x > state.max[0]) state.max[0] = x
  if (y > state.max[1]) state.max[1] = y
}

function project(x, y) {
  noteBounds(x, y)
  return [x, y]
}

/* Two-pass: first gather bounds, then recentre. To avoid storing geometry
 * twice we collect raw projected features into memory (fine for this size). */
const raw = {
  catchments: [],
  gwb: [],
  rivers: [],
  reservoirs: [],
  transitional: [],
  coastal: [],
  massesPoints: [],
}

function geomBoundsWalk(coords) {
  if (typeof coords[0] === 'number') { noteBounds(coords[0], coords[1]); return }
  for (const c of coords) geomBoundsWalk(c)
}

function walkGeometry(geom) {
  if (!geom) return
  geomBoundsWalk(geom.coordinates ?? [])
}

/* ------------------------------------------------------------------ */
/* Loaders                                                            */
/* ------------------------------------------------------------------ */

async function loadLayers() {
  for (const [key, path] of Object.entries(LAYER_FILES)) {
    if (['sites', 'connectivity', 'allMassesXlsx'].includes(key)) continue
    process.stdout.write(`  reading ${key} ... `)
    const feats = await readFeatures(path)
    raw[key] = feats
    for (const f of feats) walkGeometry(f.geometry)
    console.log(`${feats.length} features`)
  }

  process.stdout.write('  reading sites ... ')
  const siteFeats = await readFeatures(LAYER_FILES.sites)
  console.log(`${siteFeats.length} features`)

  return siteFeats
}

/* ------------------------------------------------------------------ */
/* Builders                                                           */
/* ------------------------------------------------------------------ */

function transformRing(ring, tol, ox, oy) {
  const pts = ringToPoints(ring)
  if (pts.length < 3) return null
  const simplified = simplifyLine(pts, tol)
  if (simplified.length < 3) return null
  const out = simplified.map(([x, y]) => [Math.round(x - ox), Math.round(y - oy)])
  return out
}

function transformPolygon(coords, tol, ox, oy) {
  const rings = []
  // outer ring + holes; keep holes only if they simplify to something valid
  for (let i = 0; i < coords.length; i++) {
    const r = transformRing(coords[i], i === 0 ? tol : tol, ox, oy)
    if (r && r.length >= 3) {
      if (i > 0 && ringArea(r) < 50000) continue
      rings.push(r)
    }
  }
  return rings.length ? rings : null
}

function transformGeometry(geom, tol, ox, oy) {
  if (!geom) return null
  if (geom.type === 'Polygon') {
    const rings = transformPolygon(geom.coordinates, tol, ox, oy)
    return rings ? { type: 'Polygon', coordinates: rings } : null
  }
  if (geom.type === 'MultiPolygon') {
    const polys = []
    for (const poly of geom.coordinates) {
      const rings = transformPolygon(poly, tol, ox, oy)
      if (rings) polys.push(rings)
    }
    return polys.length ? { type: 'MultiPolygon', coordinates: polys } : null
  }
  if (geom.type === 'LineString') {
    const pts = geom.coordinates.map((p) => [p[0], p[1]])
    const s = simplifyLine(pts, tol)
    if (s.length < 2) return null
    return { type: 'LineString', coordinates: s.map(([x, y]) => [Math.round(x - ox), Math.round(y - oy)]) }
  }
  if (geom.type === 'MultiLineString') {
    const lines = []
    for (const line of geom.coordinates) {
      const s = simplifyLine(line.map((p) => [p[0], p[1]]), tol)
      if (s.length >= 2) lines.push(s.map(([x, y]) => [Math.round(x - ox), Math.round(y - oy)]))
    }
    return lines.length ? { type: 'MultiLineString', coordinates: lines } : null
  }
  return null
}

function makeFeature(id, properties, geometry) {
  return { type: 'Feature', id, properties, geometry }
}

function centroidOf(geom) {
  // average of outer ring vertices, recentred later by caller
  const ring = geom.type === 'Polygon'
    ? geom.coordinates[0]
    : geom.type === 'MultiPolygon'
      ? geom.coordinates[0][0]
      : geom.coordinates
  let sx = 0, sy = 0
  for (const p of ring) { sx += p[0]; sy += p[1] }
  return [sx / ring.length, sy / ring.length]
}

function buildSimpleLayer(feats, tol, ox, oy, propsFn, filterFn) {
  const features = []
  let dropped = 0
  for (const f of feats) {
    if (filterFn && !filterFn(f)) continue
    const geometry = transformGeometry(f.geometry, tol, ox, oy)
    if (!geometry) { dropped++; continue }
    features.push(makeFeature(features.length, propsFn(f), geometry))
  }
  if (dropped) console.log(`    (dropped ${dropped} degenerate geometries)`)
  return features
}

function writeJson(name, obj) {
  return writeFile(join(OUT, name), JSON.stringify(obj), 'utf8')
}

async function sizeOf(name) {
  try {
    const s = await stat(join(OUT, name))
    return `${(s.size / 1024).toFixed(0)} KB`
  } catch { return 'n/a' }
}

/* ------------------------------------------------------------------ */
/* CSV helpers                                                        */
/* ------------------------------------------------------------------ */

function parseCsv(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.length)
  if (!lines.length) return []
  const header = lines[0].split(',')
  const rows = []
  for (let i = 1; i < lines.length; i++) {
    const cells = lines[i].split(',')
    const row = {}
    for (let j = 0; j < header.length; j++) row[header[j]] = cells[j]
    rows.push(row)
  }
  return rows
}

/* ------------------------------------------------------------------ */
/* Species dictionary                                                 */
/* ------------------------------------------------------------------ */

const SPECIES = {
  Ls: { name: 'Luciobarbus sclateri', status: 'native' },
  Ia: { name: 'Ia (see README; likely Squalius alburnoides)', status: 'native' },
  Sa: { name: 'Squalius alburnoides', status: 'native' },
  Sp: { name: 'Squalius pyrenaicus', status: 'native' },
  Pw: { name: 'Pseudochondrostoma willkommii', status: 'native' },
  Cp: { name: 'Cobitis paludica', status: 'native' },
  Il: { name: 'Iberochondrostoma lemmingii', status: 'native' },
  Io: { name: 'Iberochondrostoma oretanum', status: 'native' },
  Ah: { name: 'Anaecypris hispanica', status: 'native' },
  St: { name: 'Salmo trutta', status: 'native' },
  Lr: { name: 'Liza ramada', status: 'native' },
  Mc: { name: 'Mugil cephalus', status: 'native' },
  Ab: { name: 'Aphanius baeticus', status: 'native' },
  Gl: { name: 'Gobio lozanoi', status: 'exotic' },
  Gh: { name: 'Gambusia holbrooki', status: 'exotic' },
  Lg: { name: 'Lepomis gibbosus', status: 'exotic' },
  Aal: { name: 'Alburnus alburnus', status: 'exotic' },
  Cg: { name: 'Carassius gibelio', status: 'exotic' },
  Cc: { name: 'Cyprinus carpio', status: 'exotic' },
  Ms: { name: 'Micropterus salmoides', status: 'exotic' },
  Om: { name: 'Oncorhynchus mykiss', status: 'exotic' },
  El: { name: 'Esox lucius', status: 'exotic' },
  Am: { name: 'Ameiurus melas', status: 'exotic' },
  Tt: { name: 'Tinca tinca', status: 'exotic' },
  Aa: { name: 'Anguilla anguilla', status: 'uncertain' },
}

/* ------------------------------------------------------------------ */
/* Main                                                               */
/* ------------------------------------------------------------------ */

console.log('GuadeX data build')
await mkdir(OUT, { recursive: true })

const siteFeats = await loadLayers()

// Establish origin as the centre of the site footprint (sites are the focus).
let sxMin = Infinity, syMin = Infinity, sxMax = -Infinity, syMax = -Infinity
for (const f of siteFeats) {
  const [x, y] = f.geometry.coordinates
  sxMin = Math.min(sxMin, x); sxMax = Math.max(sxMax, x)
  syMin = Math.min(syMin, y); syMax = Math.max(syMax, y)
}
const ox = Math.round((sxMin + sxMax) / 2 / 1000) * 1000
const oy = Math.round((syMin + syMax) / 2 / 1000) * 1000
console.log(`  origin (EPSG:25830): ${ox}, ${oy}`)

// --- water-body dictionaries -------------------------------------
const massDict = {}
for (const f of raw.massesPoints) {
  const p = cleanProps(f.properties)
  if (p.ID_masa) massDict[p.ID_masa] = { name: p.Nombre_mas, system: p.ID_sistema, subzone: p.ID_subzona, zone: p.ID_zona, area: Number(p.Area_m2) || null, perimeter: Number(p.Perimetro_) || null, hasPoints: true }
}
console.log(`  water bodies with sampling points: ${Object.keys(massDict).length}`)

// --- elevations from ConnectivityUTM -----------------------------
const elevByCode = {}
try {
  const rows = parseCsv(readFileSync(LAYER_FILES.connectivity, 'utf8'))
  for (const r of rows) {
    const a = Number(r.ALTITUD)
    if (r.CODIGO && Number.isFinite(a)) elevByCode[r.CODIGO] = a
  }
} catch (e) { console.log('  warn: connectivity csv:', e.message) }
console.log(`  elevations available: ${Object.keys(elevByCode).length}`)

// --- sampling sites ----------------------------------------------
const sites = []
const siteMasaCount = {}
const seen = new Map()
for (const f of siteFeats) {
  const rawProps = cleanProps(f.properties)
  const [rx, ry] = f.geometry.coordinates
  const code = String(rawProps.CODIGO ?? '')
  if (!code) continue
  const n = (seen.get(code) ?? 0) + 1
  seen.set(code, n)
  const id = n === 1 ? code : `${code}#${n}`
  const elev = elevByCode[code]
  const masa = rawProps.ID_masa || null
  if (masa) siteMasaCount[masa] = (siteMasaCount[masa] ?? 0) + 1

  const attrs = { ...rawProps }
  delete attrs.CODIGO
  delete attrs.X
  delete attrs.Y

  sites.push({
    id,
    code,
    x: Math.round(rx - ox),
    y: Math.round(ry - oy),
    elev: Number.isFinite(elev) ? elev : null,
    subcatchment: rawProps.CODIGO_S ?? null,
    waterBody: masa,
    waterBodyName: masa && massDict[masa]?.name ? massDict[masa].name : null,
    attrs,
  })
}
console.log(`  sites: ${sites.length} (${seen.size} distinct codes)`)

// --- layers ------------------------------------------------------
console.log('  building layers ...')

const catchmentFeatures = buildSimpleLayer(
  raw.catchments, TOLERANCE.catchments, ox, oy,
  (f) => {
    const p = cleanProps(f.properties)
    const id = p.EUMASCod ?? ''
    return {
      id,
      name: p.MAS_Nombre ?? massDict[id]?.name ?? null,
      zone: p.COD_ZONA ?? null,
      subzone: p.COD_SZONA ?? null,
      areaM2: p.AREA_CUENC ? Number(p.AREA_CUENC) : null,
      perimeterM: p.PERI_CUENC ? Number(p.PERI_CUENC) : null,
      cycle: p.MASA_CICLO ?? null,
      sites: siteMasaCount[id] ?? 0,
      groundwater: false,
    }
  },
  (f) => String(cleanProps(f.properties).EUMASCod ?? '').startsWith('ES050'),
)

// Backfill the water-body dictionary from the catchment layer, then resolve any
// site water-body names that only the catchment attributes knew about.
for (const f of catchmentFeatures) {
  const p = f.properties
  if (!p.id) continue
  if (!massDict[p.id]) {
    massDict[p.id] = { name: p.name, zone: p.zone, subzone: p.subzone, area: p.areaM2, perimeter: p.perimeterM, hasPoints: false }
  } else if (!massDict[p.id].name && p.name) {
    massDict[p.id].name = p.name
  }
}
let backfilled = 0
for (const s of sites) {
  if (!s.waterBodyName && s.waterBody && massDict[s.waterBody]?.name) {
    s.waterBodyName = massDict[s.waterBody].name
    backfilled++
  }
}
console.log(`  water bodies after catchment backfill: ${Object.keys(massDict).length} (${backfilled} site names recovered)`)

const gwbFeatures = buildSimpleLayer(
  raw.gwb, TOLERANCE.gwb, ox, oy,
  (f) => {
    const p = cleanProps(f.properties)
    return {
      id: p.EUMASCod ?? p.EUDHCod ?? '',
      name: p.MAS_Nombre ?? null,
      nameEn: p.MAS_Nomb_E ?? null,
      category: p.Categoria ?? null,
      areaM2: p.AREA ? Number(p.AREA) : null,
      associated: p.Associated ?? null,
      groundwater: true,
    }
  },
  (f) => String(cleanProps(f.properties).EUDHCod ?? '').startsWith('ES050'),
)

const riverFeatures = buildSimpleLayer(
  raw.rivers, TOLERANCE.rivers, ox, oy,
  (f) => {
    const p = cleanProps(f.properties)
    return {
      id: p.EUMASCod ?? '',
      name: p.MAS_Nombre ?? null,
      category: p.Categoria ?? null,
      natural: p.Natural ?? null,
      persistence: p.Persistenc ?? null,
      meanDepth: p.meanDepth ? Number(p.meanDepth) : null,
      lengthM: p.Long ? Number(p.Long) : null,
    }
  },
  (f) => String(cleanProps(f.properties).EUDHCod ?? '').startsWith('ES050'),
)

const polygonProps = (f) => {
  const p = cleanProps(f.properties)
  return {
    id: p.EUMASCod ?? '',
    name: p.MAS_Nombre ?? null,
    category: p.Categoria ?? null,
    natural: p.Natural_1 ?? p.Natural ?? null,
    reservoir: p.reservoir ?? null,
    persistence: p.Persistenc ?? null,
    meanDepth: p.meanDepth ? Number(p.meanDepth) : null,
    areaM2: p.AREA ? Number(p.AREA) : (p.Shape_Area ? Number(p.Shape_Area) : null),
  }
}

const reservoirFeatures = buildSimpleLayer(raw.reservoirs, TOLERANCE.reservoirs, ox, oy, polygonProps, (f) => String(cleanProps(f.properties).EUDHCod ?? '').startsWith('ES050'))
const transitionalFeatures = buildSimpleLayer(raw.transitional, TOLERANCE.transitional, ox, oy, polygonProps, (f) => String(cleanProps(f.properties).EUDHCod ?? '').startsWith('ES050'))
const coastalFeatures = buildSimpleLayer(raw.coastal, TOLERANCE.coastal, ox, oy, polygonProps, (f) => String(cleanProps(f.properties).EUDHCod ?? '').startsWith('ES050'))

/* --- masked relief grid (IDW from site altitudes) ---------------- */
function buildRelief(catchmentFeatures, sites) {
  const minX = Math.round(state.min[0] - ox)
  const maxX = Math.round(state.max[0] - ox)
  const minY = Math.round(state.min[1] - oy)
  const maxY = Math.round(state.max[1] - oy)
  const cell = 1000
  const width = Math.ceil((maxX - minX) / cell) + 1
  const height = Math.ceil((maxY - minY) / cell) + 1
  const mask = new Uint8Array(width * height)

  // scanline-fill the union of catchment rings into the mask
  console.log(`  rasterising basin mask (${width}x${height}) ...`)
  for (const f of catchmentFeatures) {
    if (!f.geometry) continue
    const polys = f.geometry.type === 'Polygon' ? [f.geometry.coordinates] : f.geometry.coordinates
    for (const rings of polys) {
      for (const ring of rings) {
        for (let j = 0; j < height; j++) {
          const y = minY + j * cell
          const xs = []
          for (let k = 0; k < ring.length - 1; k++) {
            const [x1, y1] = ring[k]
            const [x2, y2] = ring[k + 1]
            if ((y1 <= y && y2 > y) || (y2 <= y && y1 > y)) {
              xs.push(x1 + ((y - y1) / (y2 - y1)) * (x2 - x1))
            }
          }
          if (xs.length < 2) continue
          xs.sort((a, b) => a - b)
          for (let s = 0; s + 1 < xs.length; s += 2) {
            const i0 = Math.max(0, Math.ceil((xs[s] - minX) / cell))
            const i1 = Math.min(width - 1, Math.floor((xs[s + 1] - minX) / cell))
            for (let i = i0; i <= i1; i++) mask[j * width + i] = 1
          }
        }
      }
    }
  }

  const pts = sites.filter((s) => s.elev != null).map((s) => [s.x, s.y, s.elev])
  console.log(`  interpolating relief from ${pts.length} altituded sites ...`)
  const data = new Int16Array(width * height).fill(-32768)
  const K = 10
  const kd = new Float64Array(K)
  const kv = new Float64Array(K)
  let maskedCells = 0
  for (let j = 0; j < height; j++) {
    const y = minY + j * cell
    for (let i = 0; i < width; i++) {
      const idx = j * width + i
      if (!mask[idx]) continue
      maskedCells++
      const x = minX + i * cell
      let count = 0
      for (let q = 0; q < K; q++) kd[q] = Infinity
      for (let p = 0; p < pts.length; p++) {
        const dx = pts[p][0] - x
        const dy = pts[p][1] - y
        const d2 = dx * dx + dy * dy
        if (count < K) {
          let m = count++
          while (m > 0 && kd[m - 1] > d2) { kd[m] = kd[m - 1]; kv[m] = kv[m - 1]; m-- }
          kd[m] = d2; kv[m] = pts[p][2]
        } else if (d2 < kd[K - 1]) {
          let m = K - 1
          while (m > 0 && kd[m - 1] > d2) { kd[m] = kd[m - 1]; kv[m] = kv[m - 1]; m-- }
          kd[m] = d2; kv[m] = pts[p][2]
        }
      }
      let sw = 0, sv = 0
      for (let q = 0; q < count; q++) {
        if (kd[q] <= 1) { sw = 1; sv = kv[q]; break }
        const w = 1 / (kd[q] * Math.sqrt(kd[q])) // 1/d^3
        sw += w; sv += w * kv[q]
      }
      data[idx] = Math.round(sv / sw)
    }
  }

  const bytes = Buffer.from(data.buffer)
  return {
    width, height, cellSize: cell, origin: [minX, minY],
    maskedCells, data: bytes.toString('base64'),
  }
}

const relief = buildRelief(catchmentFeatures, sites)
await writeJson('relief.json', relief)

const layers = {
  catchments: { file: 'catchments.json', count: catchmentFeatures.length, geometry: 'polygon' },
  gwb: { file: 'gwb.json', count: gwbFeatures.length, geometry: 'polygon' },
  rivers: { file: 'rivers.json', count: riverFeatures.length, geometry: 'line' },
  reservoirs: { file: 'reservoirs.json', count: reservoirFeatures.length, geometry: 'polygon' },
  transitional: { file: 'transitional.json', count: transitionalFeatures.length, geometry: 'polygon' },
  coastal: { file: 'coastal.json', count: coastalFeatures.length, geometry: 'polygon' },
}

await writeJson('catchments.json', { type: 'FeatureCollection', features: catchmentFeatures })
await writeJson('gwb.json', { type: 'FeatureCollection', features: gwbFeatures })
await writeJson('rivers.json', { type: 'FeatureCollection', features: riverFeatures })
await writeJson('reservoirs.json', { type: 'FeatureCollection', features: reservoirFeatures })
await writeJson('transitional.json', { type: 'FeatureCollection', features: transitionalFeatures })
await writeJson('coastal.json', { type: 'FeatureCollection', features: coastalFeatures })

// --- field metadata for the UI -----------------------------------
const elevationMin = Math.min(...sites.map((s) => s.elev).filter((v) => v != null))
const elevationMax = Math.max(...sites.map((s) => s.elev).filter((v) => v != null))

const fieldMeta = {
  code: { label: 'Site code', group: 'identity' },
  subcatchment: { label: 'Subcatchment', group: 'identity' },
  waterBody: { label: 'Water body ID', group: 'identity' },
  waterBodyName: { label: 'Water body', group: 'identity' },
  elev: { label: 'Altitude', unit: 'm', group: 'environment' },
  FECHA: { label: 'Sampling date', group: 'identity' },
  ANCHURA_m: { label: 'River width', unit: 'm', group: 'environment' },
  CODIGO_S: { label: 'Subcatchment code', group: 'identity' },
  EN_POZAS: { label: 'In pools', group: 'habitat' },
  OLOR_RESID: { label: 'Residual odour', group: 'water quality' },
  CANALIZADO: { label: 'Channelled', group: 'hydromorphology' },
  SIN_PECES: { label: 'Fishless', group: 'biota' },
  AUTOCTONAS: { label: 'Native species richness', group: 'biota' },
  ALOCTONAS: { label: 'Exotic species richness', group: 'biota' },
  TOTAL_SP: { label: 'Total species richness', group: 'biota' },
  EN: { label: 'Endangered species count', group: 'conservation' },
  VU: { label: 'Vulnerable species count', group: 'conservation' },
  LRNT: { label: 'Least concern species count', group: 'conservation' },
  ZIC: { label: 'ZIC index', group: 'conservation' },
  BCC: { label: 'BCC index', group: 'conservation' },
  VA: { label: 'VA index', group: 'conservation' },
  VA_MAX: { label: 'VA max category', group: 'conservation' },
  Causa_seco: { label: 'Drying cause', group: 'environment' },
  ID_masa: { label: 'Water body ID', group: 'identity' },
}

const manifest = {
  generated: new Date().toISOString(),
  crs: 'EPSG:25830',
  origin: [ox, oy],
  bounds: {
    minX: Math.round(state.min[0] - ox),
    minY: Math.round(state.min[1] - oy),
    maxX: Math.round(state.max[0] - ox),
    maxY: Math.round(state.max[1] - oy),
  },
  siteCount: sites.length,
  elevation: { min: elevationMin, max: elevationMax },
  relief: { file: 'relief.json', width: relief.width, height: relief.height, cellSize: relief.cellSize, origin: relief.origin },
  species: SPECIES,
  fieldMeta,
  layers,
}

await writeJson('manifest.json', manifest)
await writeJson('sites.json', { count: sites.length, sites })

// --- water bodies (333 with sampling points) as a lightweight link table
const waterBodies = Object.entries(massDict).map(([id, m]) => ({ id, ...m }))
await writeJson('waterbodies.json', { count: waterBodies.length, waterBodies })

/* --- clearly-synthetic demo results ------------------------------------
 * These exist only so the viewer's results pipeline and UI can be tested
 * end-to-end. They are deterministic but NOT real model output. */
function hashCode(s) {
  let h = 2166136261
  for (let i = 0; i < s.length; i++) { h ^= s.charCodeAt(i); h = Math.imul(h, 16777619) }
  return (h >>> 0) / 4294967295
}
const demoSteps = ['2026', '2050', '2100']
const demoSeries = {}
const demoMetrics = {}
for (const s of sites) {
  const ex = Number(s.attrs.ALOCTONAS) || 0
  const na = Number(s.attrs.AUTOCTONAS) || 0
  const elev = s.elev ?? 300
  const noise = hashCode(s.id)
  const base = Math.max(0.02, Math.min(0.8, 0.08 + 0.055 * ex - 0.00004 * elev + 0.12 * (noise - 0.5)))
  demoSeries[s.id] = demoSteps.map((_, i) => Math.round(Math.min(0.99, base * Math.pow(1.45, i)) * 1000) / 1000)
  demoMetrics[s.id] = {
    extinction_risk_2100: demoSeries[s.id][2],
    native_richness: na,
    exotic_richness: ex,
    exotic_fraction: na + ex > 0 ? Math.round((ex / (na + ex)) * 1000) / 1000 : 0,
  }
}
await writeJson('results.demo-timeseries.json', {
  name: 'DEMO — native fish extinction risk (synthetic)',
  unit: 'probability',
  description: 'SYNTHETIC demo data generated from exotic richness, elevation and a hash of the site code. Not a model result — replace with your own simulation output.',
  steps: demoSteps,
  data: demoSeries,
})
await writeJson('results.demo-metrics.json', {
  name: 'DEMO — site metrics (synthetic)',
  unit: '',
  description: 'SYNTHETIC demo data for testing the multi-metric view. Not a model result.',
  data: demoMetrics,
})

console.log('\nWritten to public/data:')
for (const f of ['manifest.json', 'sites.json', 'waterbodies.json', 'catchments.json', 'gwb.json', 'rivers.json', 'reservoirs.json', 'transitional.json', 'coastal.json', 'relief.json', 'results.demo-timeseries.json', 'results.demo-metrics.json']) {
  console.log(`  ${f.padEnd(22)} ${await sizeOf(f)}`)
}
console.log(`\nBounds (recentred m): x [${manifest.bounds.minX}, ${manifest.bounds.maxX}]  y [${manifest.bounds.minY}, ${manifest.bounds.maxY}]`)
console.log(`Elevation: ${elevationMin} - ${elevationMax} m`)
