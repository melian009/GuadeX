import { readFileSync } from 'node:fs'
import { fileURLToPath } from 'node:url'
import { dirname, join } from 'node:path'
import { Relief } from '../src/core/Relief.js'
import { normaliseResults } from '../src/data/ResultsModel.js'

const __dirname = dirname(fileURLToPath(import.meta.url))
const D = join(__dirname, '..', 'public', 'data')
const read = (f) => JSON.parse(readFileSync(join(D, f), 'utf8'))

let failures = 0
const check = (name, cond, extra = '') => {
  if (!cond) { failures++; console.log(`FAIL  ${name} ${extra}`) }
  else console.log(`ok    ${name} ${extra}`)
}

// --- relief ---------------------------------------------------------
const relief = new Relief(read('relief.json'))
const sites = read('sites.json').sites
let compared = 0
let sumAbs = 0
for (const s of sites) {
  if (s.elev == null) continue
  const h = relief.heightAt(s.x, s.y)
  if (h == null) continue
  // IDW with k=10 should stay within the local elevation envelope
  compared++
  sumAbs += Math.abs(h - s.elev)
}
check('relief covers altituded sites', compared > 500, `(${compared} sampled)`)
check('relief mean abs error < 400 m', sumAbs / compared < 400, `(MAE ${(sumAbs / compared).toFixed(1)} m)`)
check('relief min/max plausible', relief.minElev >= 0 && relief.maxElev < 3500, `(${relief.minElev}–${relief.maxElev} m)`)

// --- results: timeseries -------------------------------------------
const ts = normaliseResults(read('results.demo-timeseries.json'), 'demo')
check('timeseries mode', ts.mode === 'timeseries')
check('timeseries keys = steps', ts.keys.length === 3 && ts.keys[0].key === '2026')
check('timeseries site count', ts.values.size === sites.length, `(${ts.values.size})`)
const st = ts.stats.get('2100')
check('timeseries stats', st && st.min != null && st.max != null && !st.categorical, `(min ${st?.min} max ${st?.max})`)

// --- results: metrics ----------------------------------------------
const mt = normaliseResults(read('results.demo-metrics.json'), 'demo')
check('metrics mode', mt.mode === 'metrics')
check('metrics keys', mt.keys.length === 4, `(${mt.keys.map((k) => k.key).join(', ')})`)
check('metrics numeric', !mt.stats.get('native_richness').categorical)

// --- results: JSON array per site ---------------------------------
const arr = normaliseResults({ name: 'X', steps: ['a', 'b'], data: { '1.1.2': [1, 2] } }, 'inline')
check('inline array -> timeseries', arr.mode === 'timeseries' && arr.values.get('1.1.2').b === 2)

// --- results: CSV long --------------------------------------------
const csvLong = 'CODIGO,step,value\n1.1.2,2026,0.1\n1.1.2,2050,0.3\n1.1.3,2026,0.2\n1.1.3,2050,0.4\n'
const cl = normaliseResults(csvLong, 'long.csv')
check('csv long -> timeseries', cl.mode === 'timeseries' && cl.keys.length === 2)
check('csv long pivot', cl.values.get('1.1.2')['2050'] === 0.3)

// --- results: CSV wide with quoted field --------------------------
const csvWide = 'CODIGO,risk,population\n"1.1.2",0.31,"4200"\n1.1.3,0.22,9800\n'
const cw = normaliseResults(csvWide, 'wide.csv')
check('csv wide -> metrics', cw.mode === 'metrics' && cw.keys.length === 2)
check('csv wide values', cw.values.get('1.1.2').risk === 0.31 && cw.values.get('1.1.2').population === 4200)

// --- categorical JSON ---------------------------------------------
const cat = normaliseResults({ name: 'status', data: { a: 'high', b: 'low', c: 'high' } }, 'cat')
check('categorical detected', cat.stats.get('value').categorical === true)

console.log(failures ? `\n${failures} failure(s)` : '\nAll checks passed')
process.exit(failures ? 1 : 0)
