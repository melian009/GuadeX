import { fileURLToPath } from 'node:url'
import { dirname, join } from 'node:path'
import shapefile from 'shapefile'

const __dirname = dirname(fileURLToPath(import.meta.url))
const DATA = join(__dirname, '..', '..', 'data')

const layers = {
  points: join(DATA, 'GIS', 'capas GIS', 'capa_puntos muestreo con id masa', 'puntos muestreo con id masa.shp'),
  massesPoints: join(DATA, 'GIS', 'capas GIS', 'capa_masas con puntos muestreo', 'masas con puntos muestreo.shp'),
  cuencas: join(DATA, 'Version_02-01-2026-Masas', 'Cuencas_masas_agua_4c.shp'),
  gwb: join(DATA, 'Version_02-01-2026-Masas', 'GWB_4C.shp'),
  rivers: join(DATA, 'Version_02-01-2026-Masas', 'SW_Line_4C_.shp'),
  lago: join(DATA, 'Version_02-01-2026-Masas', 'SWB_Lago.shp'),
  trans: join(DATA, 'Version_02-01-2026-Masas', 'SWB_Transicion.shp'),
  costera: join(DATA, 'Version_02-01-2026-Masas', 'SWB_Costera.shp'),
}

async function inspectShp(label, path, maxRecords = 3) {
  try {
    const header = await shapefile.read(path)
    console.log(`\n===== ${label} =====`)
    const feats = header.features
    console.log('records:', feats.length)
    console.log('fields:', Object.keys(feats[0]?.properties ?? {}))
    for (const f of feats.slice(0, maxRecords)) {
      const p = { ...f.properties }
      for (const k of Object.keys(p)) if (typeof p[k] === 'string' && p[k].length > 60) p[k] = p[k].slice(0, 60) + '…'
      console.log(JSON.stringify(p, null, 1))
      console.log('geomType:', f.geometry?.type, 'coord0:', JSON.stringify(f.geometry?.coordinates?.[0]?.[0] ?? f.geometry?.coordinates?.[0]))
    }
  } catch (e) {
    console.log(`\n===== ${label} ERROR =====\n`, e.message)
  }
}

async function inspectShpHeaderOnly(label, path) {
  // readDbf gives fields without parsing geometry
  return new Promise((resolve) => {
    shapefile.read(path)
      .then(r => {
        const f = r.features[0]
        console.log(`\n===== ${label} =====`)
        console.log('records:', r.features.length)
        console.log('fields:', JSON.stringify(Object.keys(f?.properties ?? {})))
        resolve()
      })
      .catch(e => { console.log(label, 'ERR', e.message); resolve() })
  })
}

function inspectXlsx(label, path) {
  console.log(`\n===== XLSX ${label} =====\n  skipped (xlsx dependency removed; inspect with a spreadsheet tool)`)
}

await inspectShp('sampling points (puntos muestreo con id masa)', layers.points)
await inspectShp('masas con puntos muestreo', layers.massesPoints, 1)
await inspectShpHeaderOnly('Cuencas_masas_agua_4c', layers.cuencas)
await inspectShpHeaderOnly('GWB_4C', layers.gwb)
await inspectShpHeaderOnly('SW_Line_4C_', layers.rivers)
await inspectShpHeaderOnly('SWB_Lago', layers.lago)
await inspectShpHeaderOnly('SWB_Transicion', layers.trans)
await inspectShpHeaderOnly('SWB_Costera', layers.costera)

inspectXlsx('puntos muestreo con id', join(DATA, 'GIS', 'info atributos capas GIS', 'puntos muestreo con id.xlsx'))
inspectXlsx('info 333 masas con puntos muestreo', join(DATA, 'GIS', 'info atributos capas GIS', 'info 333 masas con puntos muestreo.xlsx'))
inspectXlsx('info todas-449 masas aguas', join(DATA, 'GIS', 'info atributos capas GIS', 'info todas-449 masas aguas.xlsx'))
