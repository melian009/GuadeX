import * as THREE from 'three'
import { SceneManager } from './core/SceneManager.js'
import { DataStore } from './data/DataStore.js'
import { normaliseResults } from './data/ResultsModel.js'
import { PolygonLayer } from './layers/PolygonLayer.js'
import { LineLayer } from './layers/LineLayer.js'
import { SitesLayer } from './layers/SitesLayer.js'
import { AppUI } from './ui/AppUI.js'
import { sceneX, sceneZ, unproject } from './lib/geo.js'
import { rampThreeColor, categoricalColor } from './lib/colors.js'

const store = new DataStore('./data/')
let sm = null
let ui = null
let relief = null
let sitesLayer = null
const layers = {}
let metricByKey = {}
let baseMetricGroups = []

const state = {
  relief: true,
  exaggeration: 14,
  heightScale: 14,
  onlyData: false,
  labels: true,
  source: 'base',
  metric: 'AUTOCTONAS',
  ramp: 'viridis',
  timeIndex: 0,
  results: null,
  selectedSiteId: null,
  selectedCatchment: -1,
}

/* ------------------------------------------------------------------ */
/* Helpers                                                            */
/* ------------------------------------------------------------------ */

function fmt(v) {
  if (v == null || v === '') return '—'
  if (typeof v === 'number') {
    if (!Number.isFinite(v)) return '—'
    if (Number.isInteger(v)) return String(v)
    const a = Math.abs(v)
    if (a !== 0 && (a < 0.001 || a >= 1e6)) return v.toExponential(2)
    return String(Math.round(v * 10000) / 10000)
  }
  return String(v)
}

function escapeHtml(s) {
  return String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]))
}

function kv(rows) {
  return `<dl class="kv">${rows.map(([k, v, cls]) => `<dt>${escapeHtml(k)}</dt><dd class="${cls ?? ''}">${v}</dd>`).join('')}</dl>`
}

function section(title, body) {
  return `<div class="info-section"><div class="info-section-title">${escapeHtml(title)}</div>${body}</div>`
}

/** Height in scene units for a geographic point, honouring the relief toggle. */
function drape(xMetres, yMetres) {
  if (!state.relief || !relief) return 0
  return relief.sceneHeightAt(xMetres, yMetres, state.exaggeration)
}

function baseValue(site, key) {
  if (key === 'elev') return site.elev
  const v = site.attrs?.[key]
  if (v == null || v === '') return null
  const n = typeof v === 'number' ? v : Number(v)
  return Number.isFinite(n) ? n : v
}

/* ------------------------------------------------------------------ */
/* Metric catalogue                                                   */
/* ------------------------------------------------------------------ */

function buildMetricCatalogue(manifest, sites) {
  const species = manifest.species ?? {}
  const presence = []
  const density = []
  const used = new Set(['elev'])
  for (const [code, info] of Object.entries(species)) {
    presence.push({ key: code, label: `${code} · ${info.name}`, unit: '', group: 'Species presence (0/1)' })
    density.push({ key: `${code}_DEN`, label: `${code} · ${info.name}`, unit: '', group: 'Species density' })
    used.add(code); used.add(`${code}_DEN`)
  }
  used.add('AUTOCTONAS'); used.add('ALOCTONAS'); used.add('TOTAL_SP')
  used.add('EN'); used.add('VU'); used.add('LRNT'); used.add('ZIC'); used.add('BCC'); used.add('VA')

  const discovered = {}
  for (const s of sites) {
    for (const [k, v] of Object.entries(s.attrs ?? {})) {
      if (used.has(k)) continue
      if (typeof v === 'number' && Number.isFinite(v)) discovered[k] = (discovered[k] ?? 0) + 1
    }
  }
  const other = Object.keys(discovered)
    .filter((k) => discovered[k] > sites.length * 0.2)
    .sort()
    .map((k) => ({ key: k, label: k, unit: '', group: 'Other numeric attributes' }))

  baseMetricGroups = [
    { label: 'Site', options: [
      { key: 'elev', label: 'Altitude', unit: 'm', group: 'Site' },
      { key: 'ANCHURA_m', label: 'River width', unit: 'm', group: 'Site' },
    ] },
    { label: 'Community', options: [
      { key: 'AUTOCTONAS', label: 'Native species richness', unit: '', group: 'Community' },
      { key: 'ALOCTONAS', label: 'Exotic species richness', unit: '', group: 'Community' },
      { key: 'TOTAL_SP', label: 'Total species richness', unit: '', group: 'Community' },
      { key: 'EN', label: 'Endangered species count', unit: '', group: 'Community' },
      { key: 'VU', label: 'Vulnerable species count', unit: '', group: 'Community' },
      { key: 'LRNT', label: 'Least-concern species count', unit: '', group: 'Community' },
    ] },
    { label: 'Indices', options: [
      { key: 'ZIC', label: 'ZIC index', unit: '', group: 'Indices' },
      { key: 'BCC', label: 'BCC index', unit: '', group: 'Indices' },
      { key: 'VA', label: 'VA index', unit: '', group: 'Indices' },
    ] },
    { label: 'Species presence', options: presence },
    { label: 'Species density', options: density },
    { label: 'Other numeric attributes', options: other },
  ]
  metricByKey = {}
  for (const g of baseMetricGroups) for (const o of g.options) metricByKey[o.key] = o
}

/** Currently active metric resolved from the source/selection. */
function metricContext() {
  if (state.source === 'results' && state.results) {
    const model = state.results
    if (model.mode === 'timeseries') {
      const step = model.keys[Math.min(state.timeIndex, model.keys.length - 1)]
      const key = step?.key
      return {
        activeKey: key,
        label: `${model.name}${key != null ? ` · ${key}` : ''}`,
        unit: model.unit ?? '',
        stats: model.globalStats ?? model.stats.get(key),
        get: (s) => model.values.get(s.id)?.[key],
        categorical: false,
      }
    }
    const chosen = model.keys.find((k) => k.key === state.metric) ?? model.keys[0]
    const key = chosen?.key
    const st = key ? model.stats.get(key) : null
    return {
      activeKey: key,
      label: chosen?.label ?? key ?? 'results',
      unit: model.unit ?? '',
      stats: st,
      get: (s) => model.values.get(s.id)?.[key],
      categorical: st?.categorical ?? state.ramp === 'categorical',
    }
  }
  const meta = metricByKey[state.metric]
  return {
    activeKey: state.metric,
    label: meta?.label ?? state.metric,
    unit: meta?.unit ?? '',
    stats: null,
    get: (s) => baseValue(s, state.metric),
    categorical: state.ramp === 'categorical',
  }
}

/* ------------------------------------------------------------------ */
/* Site style + render                                                */
/* ------------------------------------------------------------------ */

function computeStats(ctx) {
  if (ctx.stats && ctx.stats.min != null) return { min: ctx.stats.min, max: ctx.stats.max }
  let min = Infinity, max = -Infinity
  for (const s of store.sites) {
    const v = ctx.get(s)
    if (typeof v !== 'number' || !Number.isFinite(v)) continue
    if (v < min) min = v
    if (v > max) max = v
  }
  if (min === Infinity) return { min: 0, max: 1 }
  return { min, max }
}

function renderSites() {
  const ctx = metricContext()
  const numeric = typeof ctx.get(store.sites[0] ?? {}) === 'number'
  const categorical = ctx.categorical && !numeric
  const { min, max } = computeStats(ctx)
  const span = max - min || 1

  const n = store.sites.length
  const norm = new Float32Array(n)
  const raw = new Array(n)
  const hidden = new Uint8Array(n)
  for (let i = 0; i < n; i++) {
    const s = store.sites[i]
    let v = ctx.get(s)
    if (typeof v === 'string' && !categorical) {
      const num = Number(v)
      v = Number.isFinite(num) ? num : v
    }
    raw[i] = v
    if (typeof v === 'number' && Number.isFinite(v)) {
      norm[i] = (v - min) / span
    } else {
      norm[i] = NaN
    }
    if (state.onlyData && (v == null || v === '' || (typeof v === 'number' && !Number.isFinite(v)))) hidden[i] = 1
  }

  const colorFor = (i, v) => {
    if (categorical) return new THREE.Color(categoricalColor(String(raw[i] ?? '—')))
    if (v == null || Number.isNaN(v)) return new THREE.Color(0x3a4658)
    return rampThreeColor(state.ramp, v)
  }

  sitesLayer.apply(norm, colorFor, state.heightScale, hidden)
  sitesLayer.setHighlight(state.selectedSiteId ? store.sites.findIndex((s) => s.id === state.selectedSiteId) : -1)

  updateLegend(ctx, min, max, categorical, raw)
  return { ctx, min, max, raw }
}

function updateLegend(ctx, min, max, categorical, raw) {
  if (categorical) {
    const uniq = new Map()
    for (const v of raw) {
      if (v == null || v === '') continue
      const k = String(v)
      if (!uniq.has(k)) uniq.set(k, categoricalColor(k))
    }
    const categories = [...uniq.entries()].slice(0, 12).map(([label, color]) => ({ label, color }))
    ui.setLegend({ title: ctx.label, ramp: state.ramp, categorical: true, categories })
  } else {
    ui.setLegend({
      title: ctx.label,
      ramp: state.ramp,
      min: fmt(min),
      max: fmt(max),
      caption: [ctx.unit, `${store.sites.length} sites`].filter(Boolean).join(' · '),
    })
  }
}

/* ------------------------------------------------------------------ */
/* Catchment choropleth                                               */
/* ------------------------------------------------------------------ */

function updateCatchments() {
  const ctx = metricContext()
  const byWaterBody = new Map()
  for (const s of store.sites) {
    if (!s.waterBody) continue
    const v = ctx.get(s)
    if (typeof v !== 'number' || !Number.isFinite(v)) continue
    const arr = byWaterBody.get(s.waterBody) ?? []
    arr.push(v)
    byWaterBody.set(s.waterBody, arr)
  }
  const { min, max } = computeStats(ctx)
  const span = max - min || 1
  const categorical = state.ramp === 'categorical'
  const base = new THREE.Color(0x1d3347)
  const cache = new Map()

  layers.catchments.instance.setColorFunction((fi, props) => {
    if (categorical) {
      if (!props.id || !byWaterBody.has(props.id)) return base
      const vals = byWaterBody.get(props.id)
      const counts = {}
      for (const v of vals) counts[v] = (counts[v] ?? 0) + 1
      const mode = Object.entries(counts).sort((a, b) => b[1] - a[1])[0]?.[0]
      return new THREE.Color(categoricalColor(mode))
    }
    const vals = props.id ? cache.get(props.id) ?? byWaterBody.get(props.id) : null
    if (!vals || !vals.length) return base
    if (props.id) cache.set(props.id, vals)
    const mean = vals.reduce((a, b) => a + b, 0) / vals.length
    return rampThreeColor(state.ramp, (mean - min) / span)
  })
}

/* ------------------------------------------------------------------ */
/* Info panel                                                         */
/* ------------------------------------------------------------------ */

function speciesChips(site) {
  const species = store.manifest.species ?? {}
  const chips = []
  for (const [code, info] of Object.entries(species)) {
    const v = site.attrs?.[code]
    if (typeof v === 'number' && v > 0) {
      chips.push(`<span class="site-chip ${info.status}" title="${escapeHtml(info.name)} (${info.status})">${code}</span>`)
    }
  }
  return chips.length ? chips.join('') : '<span class="hint">No species recorded.</span>'
}

function densityTable(site) {
  const species = store.manifest.species ?? {}
  const rows = []
  for (const code of Object.keys(species)) {
    const v = site.attrs?.[`${code}_DEN`]
    if (typeof v === 'number' && v > 0) rows.push([`${code} · ${species[code].name}`, fmt(v), 'num'])
  }
  if (!rows.length) return ''
  return section('Species density (individuals / area)', kv(rows))
}

function infoSiteHtml(site) {
  const [utmX, utmY] = unproject(store.manifest, site.x, site.y)
  const identity = kv([
    ['Site code', escapeHtml(site.code)],
    ['Water body', escapeHtml(site.waterBodyName ?? '—')],
    ['Water body ID', escapeHtml(site.waterBody ?? '—')],
    ['Subcatchment', escapeHtml(site.subcatchment ?? '—')],
    ['Sampling date', escapeHtml(site.attrs?.FECHA ?? '—')],
    ['Altitude (m)', fmt(site.elev), 'num'],
    ['UTM easting', fmt(Math.round(utmX)), 'num'],
    ['UTM northing', fmt(Math.round(utmY)), 'num'],
  ])

  const parts = [section('Identity', identity)]

  if (state.results) {
    const model = state.results
    const rec = model.values.get(site.id)
    if (rec) {
      const rows = Object.entries(rec).map(([k, v]) => [k, escapeHtml(fmt(v)), 'num'])
      parts.push(section(`Simulation results · ${escapeHtml(model.name)}`, kv(rows)))
    } else {
      parts.push(section('Simulation results', `<div class="hint">No result value for this site.</div>`))
    }
  }

  parts.push(section('Community', kv([
    ['Native richness', fmt(site.attrs?.AUTOCTONAS), 'num'],
    ['Exotic richness', fmt(site.attrs?.ALOCTONAS), 'num'],
    ['Total richness', fmt(site.attrs?.TOTAL_SP), 'num'],
    ['ZIC', fmt(site.attrs?.ZIC), 'num'],
    ['BCC', fmt(site.attrs?.BCC), 'num'],
    ['VA', fmt(site.attrs?.VA), 'num'],
  ])))

  parts.push(section('Species present', `<div>${speciesChips(site)}</div>`))
  const dt = densityTable(site)
  if (dt) parts.push(dt)

  const skip = new Set(['FECHA', 'CODIGO_S', 'ID_masa', 'AUTOCTONAS', 'ALOCTONAS', 'TOTAL_SP', 'EN', 'VU', 'LRNT', 'ZIC', 'BCC', 'VA', 'VA_MAX', 'Causa_seco', 'ANCHURA_m'])
  for (const code of Object.keys(store.manifest.species ?? {})) { skip.add(code); skip.add(`${code}_DEN`) }
  const others = Object.entries(site.attrs ?? {})
    .filter(([k]) => !skip.has(k))
    .map(([k, v]) => [k, fmt(v)])
  if (others.length) {
    parts.push(section('All site attributes', `<dl class="kv">${others.map(([k, v]) => `<dt>${escapeHtml(k)}</dt><dd>${escapeHtml(v)}</dd>`).join('')}</dl>`))
  }

  return parts.join('')
}

function infoWaterBodyHtml(catchmentIndex) {
  const layer = layers.catchments.instance
  const props = layer.featureProps[catchmentIndex]
  const siteCount = store.sites.filter((s) => s.waterBody === props.id).length
  const parts = [section('Water body', kv([
    ['Name', escapeHtml(props.name ?? '—')],
    ['ID', escapeHtml(props.id ?? '—')],
    ['Zone', escapeHtml(props.zone ?? '—')],
    ['Sub-zone', escapeHtml(props.subzone ?? '—')],
    ['Catchment area (km²)', props.areaM2 ? fmt(Math.round(props.areaM2 / 1e6)) : '—', 'num'],
    ['Perimeter (km)', props.perimeterM ? fmt(Math.round(props.perimeterM / 100) / 10) : '—', 'num'],
    ['Sampling sites', String(siteCount), 'num'],
  ]))]
  if (siteCount) {
    const chips = store.sites
      .filter((s) => s.waterBody === props.id)
      .slice(0, 60)
      .map((s) => `<span class="site-chip" data-site="${escapeHtml(s.id)}" style="cursor:pointer">${escapeHtml(s.code)}</span>`)
      .join('')
    parts.push(section('Sites in this water body', `<div>${chips}</div>`))
  }
  return parts.join('')
}

/* ------------------------------------------------------------------ */
/* Selection                                                          */
/* ------------------------------------------------------------------ */

function selectSite(index, { focus = false } = {}) {
  if (index < 0 || index >= store.sites.length) return
  const site = store.sites[index]
  state.selectedSiteId = site.id
  if (state.selectedCatchment >= 0) {
    layers.catchments.instance.setHighlight(-1)
    state.selectedCatchment = -1
  }
  sitesLayer.setHighlight(index)
  ui.setSearchValue(site.id)
  ui.showInfo({ title: site.code, subtitle: site.waterBodyName ?? 'Sampling site', html: infoSiteHtml(site) })
  if (focus) {
    const y = sitesLayer.baseY[index] ?? 0
    sm.focus(new THREE.Vector3(sceneX(site.x), y + 4, sceneZ(site.y)), 26)
  }
}

function selectCatchment(featureIndex) {
  const layer = layers.catchments.instance
  if (featureIndex < 0) return
  state.selectedCatchment = featureIndex
  state.selectedSiteId = null
  sitesLayer.setHighlight(-1)
  layer.setHighlight(featureIndex, new THREE.Color(0xffffff))
  const props = layer.featureProps[featureIndex]
  const anySite = store.sites.find((s) => s.waterBody === props.id)
  const sub = props.name ?? 'Water body'
  ui.showInfo({ title: props.name ?? props.id ?? 'Water body', subtitle: sub, html: infoWaterBodyHtml(featureIndex) })
  if (anySite) ui.setSearchValue(anySite.id)
}

/* ------------------------------------------------------------------ */
/* Results loading                                                    */
/* ------------------------------------------------------------------ */

function rebuildMetricOptions() {
  if (state.source === 'results' && state.results) {
    const model = state.results
    if (model.mode === 'timeseries') {
      ui.setMetricOptions([{ label: 'Simulation results', options: [{ key: '__series__', label: `${model.name} (time series)` }] }], '__series__')
    } else {
      ui.setMetricOptions([{ label: 'Simulation results', options: model.keys.map((k) => ({ key: k.key, label: k.label })) }], state.metric)
    }
  } else {
    ui.setMetricOptions(baseMetricGroups, state.metric)
  }
}

/** Keep the selected metric valid for the active source. */
function ensureMetricForSource() {
  if (state.source === 'results' && state.results) {
    if (state.results.mode === 'timeseries') state.metric = '__series__'
    else if (!state.results.keys.some((k) => k.key === state.metric)) state.metric = state.results.keys[0]?.key ?? 'value'
  } else if (!metricByKey[state.metric]) {
    state.metric = 'AUTOCTONAS'
  }
}

function loadResultsPayload(payload, name) {
  let model
  try {
    model = normaliseResults(payload, name)
  } catch (e) {
    ui.setResultsStatus(`Could not parse results: ${e.message}`, 'err')
    return
  }
  if (!model.values.size) {
    ui.setResultsStatus('No sites found in the results file.', 'err')
    return
  }
  state.results = model
  state.source = 'results'
  ui.el.selSource.value = 'results'
  state.timeIndex = 0
  state.metric = model.mode === 'timeseries' ? '__series__' : (model.keys[0]?.key ?? 'value')
  ensureMetricForSource()
  rebuildMetricOptions()
  ui.setTime(model.mode === 'timeseries' ? model.keys.map((k) => k.label) : null, 0)
  ui.setResultsStatus(`Loaded “${model.name}” — ${model.values.size} sites, ${model.keys.length} ${model.mode === 'timeseries' ? 'steps' : 'metrics'}.`, 'ok')
  renderAll()
}

async function loadResultsFromUrl(url) {
  const res = await fetch(url)
  if (!res.ok) throw new Error(`${res.status} ${res.statusText}`)
  const text = await res.text()
  loadResultsPayload(text, url.split('/').pop())
}

function clearResults() {
  state.results = null
  state.source = 'base'
  ui.el.selSource.value = 'base'
  state.timeIndex = 0
  ensureMetricForSource()
  rebuildMetricOptions()
  ui.setTime(null, 0)
  ui.setResultsStatus('No results loaded — showing site attributes.')
  renderAll()
}

async function handleFile(file) {
  const text = await file.text()
  loadResultsPayload(text, file.name)
}

/* ------------------------------------------------------------------ */
/* Render orchestration                                               */
/* ------------------------------------------------------------------ */

let lastCtx = null

function renderAll() {
  if (!sitesLayer) return
  const info = renderSites()
  lastCtx = info.ctx
  updateCatchments()
  updateStats(info)
}

function updateStats(info) {
  ui.setStats([
    { label: 'sites', value: store.sites.length },
    { label: 'water bodies', value: layers.catchments?.instance.features.length ?? 0 },
    { label: 'metric', value: info.ctx.label },
    { label: 'range', value: `${fmt(info.min)} – ${fmt(info.max)}` },
    { label: state.results ? 'results' : 'source', value: state.results ? state.results.name : 'attributes' },
  ])
}

const labelEl = typeof document !== 'undefined' ? document.getElementById('site-label') : null
const _labelVec = new THREE.Vector3()

function updateLabel() {
  if (!labelEl || !sitesLayer) return
  const idx = state.selectedSiteId ? store.sites.findIndex((s) => s.id === state.selectedSiteId) : -1
  if (!state.labels || idx < 0) { labelEl.classList.add('hidden'); return }
  const s = store.sites[idx]
  _labelVec.set(sceneX(s.x), (sitesLayer.baseY[idx] ?? 0) + sitesLayer.heights[idx] + 2.2, sceneZ(s.y))
  const p = sm.projectToScreen(_labelVec)
  if (!p.visible) { labelEl.classList.add('hidden'); return }
  labelEl.textContent = s.code
  labelEl.style.left = `${p.x}px`
  labelEl.style.top = `${p.y}px`
  labelEl.classList.remove('hidden')
}

/* ------------------------------------------------------------------ */
/* Scene construction                                                 */
/* ------------------------------------------------------------------ */

async function buildScene(data) {
  const manifest = data.manifest
  relief = data.relief

  const reliefMesh = relief.buildMesh()
  relief.applyExaggeration(state.exaggeration)
  reliefMesh.visible = state.relief
  sm.add(reliefMesh)

  // catchments (base choropleth layer)
  layers.catchments = {
    instance: new PolygonLayer({
      name: 'catchments',
      features: data.features.catchments,
      color: 0x24425c,
      opacity: 0.95,
      heightAt: drape,
      outlineColor: 0x74a8d6,
      outlineOpacity: 0.35,
      renderOrder: 0,
    }),
    visible: true,
  }
  sm.add(layers.catchments.instance.mesh)
  sm.add(layers.catchments.instance.outlineMesh)
  sm.add(layers.catchments.instance.pickTarget)

  layers.gwb = {
    instance: new PolygonLayer({
      name: 'gwb',
      features: data.features.gwb,
      color: 0x5a4a9e,
      opacity: 0.45,
      heightAt: drape,
      outlineColor: 0x8f7fe0,
      outlineOpacity: 0.4,
      renderOrder: 0.1,
    }),
    visible: false,
  }
  sm.add(layers.gwb.instance.mesh)
  sm.add(layers.gwb.instance.outlineMesh)
  sm.add(layers.gwb.instance.pickTarget)

  layers.rivers = {
    instance: new LineLayer({
      name: 'rivers',
      features: data.features.rivers,
      color: 0x6cc6ff,
      width: 1.7,
      opacity: 0.9,
      heightAt: drape,
      renderOrder: 2,
    }),
    visible: true,
  }
  sm.add(layers.rivers.instance.line)

  layers.reservoirs = {
    instance: new PolygonLayer({
      name: 'reservoirs',
      features: data.features.reservoirs,
      color: 0x1f7fa8,
      opacity: 0.9,
      heightAt: drape,
      outlineColor: 0x5fd0ef,
      outlineOpacity: 0.6,
      renderOrder: 3,
    }),
    visible: true,
  }
  sm.add(layers.reservoirs.instance.mesh)
  sm.add(layers.reservoirs.instance.outlineMesh)
  sm.add(layers.reservoirs.instance.pickTarget)

  layers.transitional = {
    instance: new PolygonLayer({
      name: 'transitional',
      features: data.features.transitional,
      color: 0x1f8f78,
      opacity: 0.8,
      outlineColor: 0x53d6b4,
      outlineOpacity: 0.6,
      renderOrder: 3,
    }),
    visible: false,
  }
  sm.add(layers.transitional.instance.mesh)
  sm.add(layers.transitional.instance.outlineMesh)
  sm.add(layers.transitional.instance.pickTarget)

  layers.coastal = {
    instance: new PolygonLayer({
      name: 'coastal',
      features: data.features.coastal,
      color: 0x2f5fa8,
      opacity: 0.75,
      outlineColor: 0x6f9be0,
      outlineOpacity: 0.5,
      renderOrder: 3,
    }),
    visible: false,
  }
  sm.add(layers.coastal.instance.mesh)
  sm.add(layers.coastal.instance.outlineMesh)
  sm.add(layers.coastal.instance.pickTarget)

  // sites
  sitesLayer = new SitesLayer({ sites: store.sites })
  sitesLayer.setBaseFromRelief(relief, state.exaggeration)
  sm.add(sitesLayer.columns)
  sm.add(sitesLayer.caps)

  registerPickables()
  setLayerVisibility()
  sm.fitBounds(manifest.bounds, { padding: 1.25 })

  buildMetricCatalogue(manifest, store.sites)
  rebuildMetricOptions()
  ui.setSiteCodes(store.sites.map((s) => s.id))
  ui.setTime(null, 0)
  renderAll()
}

function registerPickables() {
  const catchL = layers.catchments.instance
  sm.registerPickable(catchL.pickTarget, {
    hitKey: (hit) => `catch:${catchL.featureIndexAt(hit)}`,
    onHover: (hit, ev) => {
      if (!hit) { ui.hideTooltip(); return }
      const fi = catchL.featureIndexAt(hit)
      if (fi < 0) { ui.hideTooltip(); return }
      const p = catchL.featureProps[fi]
      ui.showTooltip(`<b>${escapeHtml(p.name ?? p.id ?? 'Water body')}</b><br>${p.sites ?? 0} sampling sites<br><span style="color:#93a4bb">${escapeHtml(p.id ?? '')}</span>`, ev.clientX, ev.clientY)
    },
    onMove: (hit, ev) => { if (hit) ui.showTooltip(ui.el.tooltip.innerHTML, ev.clientX, ev.clientY) },
    onClick: (hit) => {
      const fi = catchL.featureIndexAt(hit)
      if (fi >= 0) selectCatchment(fi)
    },
  })

  const gwbL = layers.gwb.instance
  sm.registerPickable(gwbL.pickTarget, {
    hitKey: (hit) => `gwb:${gwbL.featureIndexAt(hit)}`,
    onHover: (hit, ev) => {
      if (!hit) { ui.hideTooltip(); return }
      const fi = gwbL.featureIndexAt(hit)
      if (fi < 0) { ui.hideTooltip(); return }
      const p = gwbL.featureProps[fi]
      ui.showTooltip(`<b>${escapeHtml(p.name ?? p.id ?? 'Groundwater body')}</b><br><span style="color:#93a4bb">${escapeHtml(p.category ?? '')}</span>`, ev.clientX, ev.clientY)
    },
  })

  const resL = layers.reservoirs.instance
  sm.registerPickable(resL.pickTarget, {
    hitKey: (hit) => `res:${resL.featureIndexAt(hit)}`,
    onHover: (hit, ev) => {
      if (!hit) { ui.hideTooltip(); return }
      const fi = resL.featureIndexAt(hit)
      if (fi < 0) { ui.hideTooltip(); return }
      const p = resL.featureProps[fi]
      ui.showTooltip(`<b>${escapeHtml(p.name ?? p.id ?? 'Water body')}</b><br><span style="color:#93a4bb">${escapeHtml(p.category ?? 'Reservoir / lake')}</span>`, ev.clientX, ev.clientY)
    },
  })

  sm.registerPickable(sitesLayer.columns, {
    hitKey: (hit) => `site:${sitesLayer.siteIndexAt(hit)}`,
    onHover: (hit, ev) => {
      if (!hit) { ui.hideTooltip(); return }
      const idx = sitesLayer.siteIndexAt(hit)
      if (idx < 0 || sitesLayer.hidden[idx]) { ui.hideTooltip(); return }
      const s = store.sites[idx]
      const ctx = lastCtx ?? metricContext()
      const v = ctx.get(s)
      ui.showTooltip(`<b>${escapeHtml(s.code)}</b><br>${escapeHtml(s.waterBodyName ?? '')}<br><span style="color:#93a4bb">${escapeHtml(ctx.label)}:</span> ${escapeHtml(fmt(v))}`, ev.clientX, ev.clientY)
    },
    onMove: (hit, ev) => { if (hit) ui.showTooltip(ui.el.tooltip.innerHTML, ev.clientX, ev.clientY) },
    onClick: (hit) => {
      const idx = sitesLayer.siteIndexAt(hit)
      if (idx >= 0 && !sitesLayer.hidden[idx]) selectSite(idx)
    },
  })
}

function setLayerVisibility() {
  for (const id of Object.keys(layers)) {
    const l = layers[id]
    if (!l) continue
    l.instance.setVisible(l.visible)
  }
}

function refreshHeights() {
  for (const id of ['catchments', 'gwb', 'reservoirs', 'transitional', 'coastal']) {
    layers[id]?.instance.refreshHeights?.()
  }
  layers.rivers?.instance.refreshHeights?.()
  relief?.setVisible(state.relief)
  if (sitesLayer) {
    sitesLayer.setBaseFromRelief(relief, state.exaggeration)
    renderAll()
  }
}

/** Coalesce exaggeration slider input to one rebuild per animation frame. */
let exaggRaf = 0
function scheduleExaggeration() {
  if (exaggRaf) return
  exaggRaf = requestAnimationFrame(() => {
    exaggRaf = 0
    relief?.applyExaggeration(state.exaggeration)
    refreshHeights()
  })
}

/* ------------------------------------------------------------------ */
/* Handlers + boot                                                    */
/* ------------------------------------------------------------------ */

const handlers = {
  onLayerToggle: (id, visible) => {
    if (!layers[id]) return
    layers[id].visible = visible
    layers[id].instance.setVisible(visible)
    if (id === 'catchments') renderAll()
  },
  onReliefToggle: (v) => {
    state.relief = v
    refreshHeights()
  },
  onExaggeration: (v) => {
    state.exaggeration = v
    scheduleExaggeration()
  },
  onHeightScale: (v) => { state.heightScale = v; renderAll() },
  onOnlyData: (v) => { state.onlyData = v; renderAll() },
  onLabels: (v) => { state.labels = v },
  onSourceChange: (v) => {
    state.source = v
    if (v === 'results' && !state.results) {
      ui.setResultsStatus('No results loaded yet — load a file to use this source.', 'err')
    }
    ensureMetricForSource()
    rebuildMetricOptions()
    renderAll()
  },
  onMetricChange: (key) => { state.metric = key; renderAll() },
  onRampChange: (v) => { state.ramp = v; renderAll() },
  onTimeChange: (i) => {
    state.timeIndex = i
    if (state.results?.mode === 'timeseries') {
      const label = state.results.keys[i]?.label
      ui.el.outTime.textContent = label ?? '—'
    }
    renderAll()
  },
  onLoadFile: (file) => { handleFile(file).catch((e) => ui.setResultsStatus(`Load failed: ${e.message}`, 'err')) },
  onLoadDemo: (kind) => {
    const url = kind === 'timeseries' ? './data/results.demo-timeseries.json' : './data/results.demo-metrics.json'
    loadResultsFromUrl(url).catch((e) => ui.setResultsStatus(`Demo load failed: ${e.message}`, 'err'))
  },
  onClearResults: clearResults,
  onSearch: (value) => {
    const v = String(value || '').trim()
    if (!v) return
    const i = store.sites.findIndex((s) => s.id === v || s.code === v)
    if (i >= 0) selectSite(i, { focus: true })
  },
  onResetView: () => sm.resetView(),
  onInfoClose: () => {},
}

async function boot() {
  sm = new SceneManager(document.getElementById('scene'))
  ui = new AppUI(handlers)

  // clicking a site chip inside the water-body panel selects it
  ui.el.infoBody.addEventListener('click', (e) => {
    const chip = e.target.closest('[data-site]')
    if (!chip) return
    const i = store.sites.findIndex((s) => s.id === chip.dataset.site)
    if (i >= 0) selectSite(i, { focus: true })
  })

  try {
    const data = await store.load((msg) => ui.setLoadingText(msg))
    await buildScene(data)
    sm.onResize((w, h) => layers.rivers?.instance.setResolution(w, h))
    sm.onFrame(updateLabel)
    ui.fadeLoading()
  } catch (e) {
    ui.setLoadingText(`Failed to load data: ${e.message}`)
    console.error(e)
    return
  }

  window.GuadeX = {
    version: '0.1.0',
    state,
    viewer: sm,
    store,
    setResults: (obj) => loadResultsPayload(obj, 'API results'),
    loadResults: (url) => loadResultsFromUrl(url),
    clearResults,
    selectSite: (id) => {
      const i = store.sites.findIndex((s) => s.id === id || s.code === id)
      if (i >= 0) selectSite(i, { focus: true })
      return i >= 0
    },
    setMetric: (key, source = 'base') => {
      state.source = source
      state.metric = key
      ui.el.selSource.value = source
      ensureMetricForSource()
      rebuildMetricOptions()
      if (state.source === 'results' && state.results?.mode === 'timeseries') state.timeIndex = 0
      renderAll()
    },
    setRamp: (name) => { state.ramp = name; ui.el.selRamp.value = name; renderAll() },
    setTimeIndex: (i) => { state.timeIndex = i; ui.el.rangeTime.value = String(i); renderAll() },
  }

  // shareable deep links:
  //   ?results=<url>&metric=<key>&site=<code>&ramp=<name>
  const params = new URLSearchParams(location.search)
  try {
    if (params.get('results')) await loadResultsFromUrl(params.get('results'))
    if (params.get('ramp')) window.GuadeX.setRamp(params.get('ramp'))
    if (params.get('metric')) window.GuadeX.setMetric(params.get('metric'), params.get('source') ?? (params.get('results') ? 'results' : 'base'))
    if (params.get('site')) window.GuadeX.selectSite(params.get('site'))
  } catch (e) {
    ui.setResultsStatus(`Deep link failed: ${e.message}`, 'err')
  }

  window.dispatchEvent(new CustomEvent('guadex:ready'))
}

boot()
