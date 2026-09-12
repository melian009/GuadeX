/**
 * Normalises heterogeneous simulation outputs into one uniform model:
 *
 *   {
 *     name, unit, description, source,
 *     keys:   [{ key, label, unit }],
 *     values: Map<siteId, Record<key, number|string>>,
 *     stats:  Map<key, { min, max, mean, count, unique, categorical }>,
 *     raw:    the original payload (shown verbatim in the site panel),
 *   }
 */

const ID_CANDIDATES = ['codigo', 'id', 'site', 'site_id', 'siteid', 'code', 'station']
const STEP_CANDIDATES = ['step', 'time', 'year', 'scenario', 'date', 't', 'period']

export function parseCsv(text) {
  const rows = []
  let row = []
  let field = ''
  let quoted = false
  for (let i = 0; i < text.length; i++) {
    const c = text[i]
    if (quoted) {
      if (c === '"') {
        if (text[i + 1] === '"') { field += '"'; i++ } else quoted = false
      } else field += c
    } else if (c === '"') {
      quoted = true
    } else if (c === ',') {
      row.push(field); field = ''
    } else if (c === '\n') {
      row.push(field); rows.push(row); row = []; field = ''
    } else if (c === '\r') {
      // ignore
    } else field += c
  }
  if (field.length || row.length) { row.push(field); rows.push(row) }
  if (!rows.length) return []
  const header = rows[0].map((h) => h.trim())
  return rows.slice(1)
    .filter((r) => r.some((c) => c !== ''))
    .map((r) => {
      const o = {}
      header.forEach((h, i) => { o[h] = (r[i] ?? '').trim() })
      return o
    })
}

function toNumber(v) {
  if (v == null || v === '') return null
  if (typeof v === 'number') return Number.isFinite(v) ? v : null
  const n = Number(String(v).replace(',', '.'))
  return Number.isFinite(n) ? n : null
}

function detectIdColumn(headers) {
  const lower = headers.map((h) => h.toLowerCase())
  for (const cand of ID_CANDIDATES) {
    const i = lower.indexOf(cand)
    if (i >= 0) return headers[i]
  }
  return headers[0]
}

function detectStepColumn(headers) {
  const lower = headers.map((h) => h.toLowerCase())
  for (const cand of STEP_CANDIDATES) {
    const i = lower.indexOf(cand)
    if (i >= 0) return headers[i]
  }
  return null
}

function finiteStats(model) {
  for (const { key } of model.keys) {
    let min = Infinity, max = -Infinity, sum = 0, count = 0
    let sawNumber = false
    let sawOther = false
    const unique = new Set()
    for (const rec of model.values.values()) {
      const v = rec[key]
      if (v == null) continue
      if (typeof v === 'number') {
        sawNumber = true
        if (v < min) min = v
        if (v > max) max = v
        sum += v
        count++
        if (unique.size <= 24) unique.add(v)
      } else {
        sawOther = true
        unique.add(v)
      }
    }
    model.stats.set(key, {
      min: count ? min : null,
      max: count ? max : null,
      mean: count ? sum / count : null,
      count,
      unique: [...unique],
      categorical: !sawNumber && sawOther,
    })
  }
  return model
}

function fromJson(payload, source) {
  const name = payload.name ?? payload.title ?? source ?? 'Results'
  const unit = payload.unit ?? ''
  const description = payload.description ?? ''
  const data = payload.data ?? payload.values ?? payload
  if (data == null || typeof data !== 'object') throw new Error('No data object found')

  const keys = []
  const values = new Map()
  const addKey = (k) => { if (!keys.includes(k)) keys.push(k) }
  let sawArray = false

  // reuse declared steps as key labels where possible
  const steps = Array.isArray(payload.steps) ? payload.steps.map(String) : null

  for (const [siteId, value] of Object.entries(data)) {
    if (siteId === 'meta' || siteId === 'name' || siteId === 'unit') continue
    const rec = {}
    if (value == null) {
      values.set(siteId, rec)
      continue
    }
    if (typeof value === 'number' || typeof value === 'string') {
      const k = payload.valueKey ?? 'value'
      addKey(k)
      rec[k] = typeof value === 'string' ? value : value
    } else if (Array.isArray(value)) {
      sawArray = true
      value.forEach((v, i) => {
        const k = steps && steps[i] != null ? steps[i] : String(i)
        addKey(k)
        rec[k] = v
      })
    } else if (typeof value === 'object') {
      for (const [k, v] of Object.entries(value)) {
        addKey(k)
        rec[k] = v
      }
    }
    values.set(siteId, rec)
  }
  return finiteStats({ name, unit, description, source, mode: sawArray ? 'timeseries' : 'metrics', steps, keys: keys.map((k) => ({ key: k, label: k, unit: '' })), values, stats: new Map(), raw: payload })
}

function fromCsv(text, source) {
  const rows = parseCsv(text)
  if (!rows.length) throw new Error('Empty CSV')
  const headers = Object.keys(rows[0])
  const idCol = detectIdColumn(headers)
  const stepCol = detectStepColumn(headers)
  const numericCols = headers.filter((h) => h !== idCol && (stepCol == null || h !== stepCol) && rows.some((r) => toNumber(r[h]) != null))

  const values = new Map()
  const keys = []
  const addKey = (k) => { if (!keys.includes(k)) keys.push(k) }

  if (stepCol) {
    for (const r of rows) {
      const id = String(r[idCol])
      const step = String(r[stepCol])
      const rec = values.get(id) ?? {}
      for (const col of numericCols) {
        const v = toNumber(r[col])
        if (v == null) continue
        const key = numericCols.length === 1 ? step : `${col} @ ${step}`
        addKey(key)
        rec[key] = v
      }
      values.set(id, rec)
    }
  } else {
    for (const r of rows) {
      const id = String(r[idCol])
      const rec = values.get(id) ?? {}
      for (const col of numericCols) {
        const v = toNumber(r[col])
        if (v != null) { rec[col] = v; addKey(col) }
      }
      values.set(id, rec)
    }
  }
  return finiteStats({ name: source ?? 'CSV results', unit: '', description: '', source, mode: stepCol ? 'timeseries' : 'metrics', keys: keys.map((k) => ({ key: k, label: k, unit: '' })), values, stats: new Map(), raw: { rows: rows.length, headers } })
}

export function normaliseResults(input, source) {
  if (typeof input === 'string') {
    const trimmed = input.trim()
    if (trimmed.startsWith('{') || trimmed.startsWith('[')) {
      return fromJson(JSON.parse(trimmed), source)
    }
    return fromCsv(trimmed, source)
  }
  return fromJson(input, source)
}
