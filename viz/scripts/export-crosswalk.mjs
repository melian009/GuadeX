/*
 * GuadeX site -> reporting-level crosswalk.
 * -----------------------------------------
 * Exports a small, dependency-free CSV that maps every sampling site to its
 * four-level reporting hierarchy:
 *
 *   sampling point (CODIGO)
 *     -> sub-basin  (CODIGO_S, the model subcatchment)
 *     -> water body (ID_masa / EUMASCod, with name, zone, subzone)
 *     -> whole basin (ES050, implied by the presence of a water body)
 *
 * The source of truth is the GIS archive:
 *   - sites:      data/GIS/capas GIS/capa_puntos muestreo con id masa/*.shp
 *   - masses:     data/GIS/capas GIS/capa_masas con puntos muestreo/*.shp
 *   - catchments: data/Version_02-01-2026-Masas/Cuencas_masas_agua_4c.shp
 *
 * Only attribute tables are read (no geometry), and the shapefile text repair
 * is shared with `build-data.mjs` through `gis-utils.mjs` so names stay in
 * lockstep with the viewer's `sites.json`.
 *
 * Run with:  npm run crosswalk
 */
import { fileURLToPath } from 'node:url'
import { dirname, join } from 'node:path'
import { writeFile } from 'node:fs/promises'
import { fixText, readFeatures, csvCell } from './gis-utils.mjs'

const __dirname = dirname(fileURLToPath(import.meta.url))
const ROOT = join(__dirname, '..')
const DATA = join(ROOT, '..', 'data')
const OUT = join(DATA, 'site_waterbody_crosswalk.csv')

const SITES_SHP = join(DATA, 'GIS', 'capas GIS', 'capa_puntos muestreo con id masa', 'puntos muestreo con id masa.shp')
const MASSES_SHP = join(DATA, 'GIS', 'capas GIS', 'capa_masas con puntos muestreo', 'masas con puntos muestreo.shp')
const CATCHMENTS_SHP = join(DATA, 'Version_02-01-2026-Masas', 'Cuencas_masas_agua_4c.shp')

function text(value) {
  if (value === null || value === undefined) return ''
  return fixText(String(value)).trim()
}

console.log('GuadeX site crosswalk export')

// --- water-body dictionary ---------------------------------------------
const mass = {}
const massFeatures = await readFeatures(MASSES_SHP)
for (const f of massFeatures) {
  const p = f.properties
  const id = text(p.ID_masa ?? p.EUMASCod)
  if (!id) continue
  mass[id] = {
    name: text(p.Nombre_mas ?? p.MAS_Nombre),
    subzone: text(p.ID_subzona),
    zone: text(p.ID_zona),
  }
}
console.log(`  water bodies with sampling points: ${Object.keys(mass).length}`)

// --- backfill names/hierarchy from the catchment layer ------------------
const catchmentFeatures = await readFeatures(CATCHMENTS_SHP)
let backfilled = 0
for (const f of catchmentFeatures) {
  const p = f.properties
  const id = text(p.EUMASCod)
  if (!id || !String(id).startsWith('ES050')) continue
  const name = text(p.MAS_Nombre)
  if (!mass[id]) {
    mass[id] = { name, subzone: text(p.COD_SZONA), zone: text(p.COD_ZONA) }
  } else if (!mass[id].name && name) {
    mass[id].name = name
    backfilled++
  }
}
console.log(`  water bodies after catchment backfill: ${Object.keys(mass).length} (${backfilled} names recovered)`)

// --- site rows ----------------------------------------------------------
const siteFeatures = await readFeatures(SITES_SHP)
const rows = []
const seen = new Set()
for (const f of siteFeatures) {
  const p = f.properties
  const code = text(p.CODIGO)
  if (!code) continue
  const subcatchment = text(p.CODIGO_S)
  const waterBody = text(p.ID_masa)
  const meta = mass[waterBody] ?? {}
  const key = [code, subcatchment, waterBody].join('|')
  if (seen.has(key)) continue
  seen.add(key)
  rows.push({
    CODIGO: code,
    CODIGO_S: subcatchment,
    ID_masa: waterBody,
    water_body_name: meta.name ?? '',
    subzone: meta.subzone ?? '',
    zone: meta.zone ?? '',
  })
}
console.log(`  sites: ${siteFeatures.length} features, ${rows.length} unique rows`)

const header = ['CODIGO', 'CODIGO_S', 'ID_masa', 'water_body_name', 'subzone', 'zone']
const lines = [header.join(',')]
for (const r of rows) lines.push(header.map((h) => csvCell(r[h])).join(','))
await writeFile(OUT, lines.join('\n') + '\n', 'utf8')
console.log(`  wrote ${OUT} (${rows.length} rows)`)
