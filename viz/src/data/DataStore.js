import { Relief } from '../core/Relief.js'

async function fetchJson(url) {
  const res = await fetch(url)
  if (!res.ok) throw new Error(`${res.status} ${res.statusText} — ${url}`)
  return res.json()
}

/**
 * Loads every pre-built dataset and exposes convenient lookups.
 */
export class DataStore {
  constructor(base = './data/') {
    this.base = base
    this.manifest = null
    this.sites = []
    this.siteById = new Map()
    this.attributes = null
  }

  async load(onProgress = () => {}) {
    onProgress('Loading manifest…')
    this.manifest = await fetchJson(this.base + 'manifest.json')

    const names = ['catchments', 'gwb', 'rivers', 'reservoirs', 'transitional', 'coastal']
    onProgress('Loading GIS layers…')
    const collections = await Promise.all(
      names.map((n) => fetchJson(this.base + this.manifest.layers[n].file)),
    )
    const byName = {}
    names.forEach((n, i) => { byName[n] = collections[i].features })

    onProgress('Loading sampling sites…')
    const sitesJson = await fetchJson(this.base + 'sites.json')
    this.sites = sitesJson.sites
    for (const s of this.sites) this.siteById.set(s.id, s)

    onProgress('Loading relief…')
    const reliefJson = await fetchJson(this.base + 'relief.json')
    const relief = new Relief(reliefJson)

    return {
      manifest: this.manifest,
      relief,
      features: byName,
      sites: this.sites,
    }
  }

  getSite(id) {
    return this.siteById.get(id) ?? null
  }

  siteIds() {
    return this.sites.map((s) => s.id)
  }
}
