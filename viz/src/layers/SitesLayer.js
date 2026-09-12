import * as THREE from 'three'
import { sceneX, sceneZ, SCALE } from '../lib/geo.js'

/**
 * All 1,037 sampling sites as GPU-instanced columns with per-instance colour.
 * Column height and colour encode a selectable variable.
 */
export class SitesLayer {
  constructor({ sites, radius = 0.7, color = 0x6fe3cd }) {
    this.sites = sites
    this.count = sites.length
    this.defaultColor = new THREE.Color(color)
    this.heights = new Float32Array(this.count).fill(3)
    this.baseY = new Float32Array(this.count)
    this.colors = new Float32Array(this.count * 3)
    this.highlightIndex = -1
    this.highlightColor = new THREE.Color(0xffffff)
    this.hidden = new Uint8Array(this.count)

    // scratch objects reused every update to avoid per-frame allocation churn
    this._c = new THREE.Color()
    this._v = new THREE.Vector3()
    this._s = new THREE.Vector3()
    this._q = new THREE.Quaternion()
    this._m = new THREE.Matrix4()

    const columnGeo = new THREE.CylinderGeometry(radius, radius, 1, 7, 1, false)
    columnGeo.translate(0, 0.5, 0)
    const columnMat = new THREE.MeshStandardMaterial({ roughness: 0.55, metalness: 0.05 })
    this.columns = new THREE.InstancedMesh(columnGeo, columnMat, this.count)
    this.columns.name = 'sites-columns'
    this.columns.instanceMatrix.setUsage(THREE.DynamicDrawUsage)
    this.columns.frustumCulled = false

    const capGeo = new THREE.SphereGeometry(radius * 1.5, 8, 6)
    const capMat = new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0.1 })
    this.caps = new THREE.InstancedMesh(capGeo, capMat, this.count)
    this.caps.name = 'sites-caps'
    this.caps.instanceMatrix.setUsage(THREE.DynamicDrawUsage)
    this.caps.frustumCulled = false

    for (let i = 0; i < this.count; i++) {
      this.setColor(i, this.defaultColor.r, this.defaultColor.g, this.defaultColor.b)
    }
    this._writeMatrices()
  }

  get pickTarget() { return this.columns }

  setColor(i, r, g, b) {
    this.colors[i * 3] = r
    this.colors[i * 3 + 1] = g
    this.colors[i * 3 + 2] = b
    this._c.setRGB(r, g, b)
    this.columns.setColorAt(i, this._c)
    this.caps.setColorAt(i, this._c)
  }

  _writeMatrices() {
    const m = this._m
    const q = this._q
    const v = this._v
    const s = this._s
    for (let i = 0; i < this.count; i++) {
      const site = this.sites[i]
      const x = sceneX(site.x)
      const z = sceneZ(site.y)
      const y = this.baseY[i]
      const h = Math.max(0.05, this.heights[i])
      const hidden = this.hidden[i] === 1
      v.set(x, y, z)
      if (hidden) s.set(0.0001, 0.0001, 0.0001); else s.set(1, h, 1)
      m.compose(v, q, s)
      this.columns.setMatrixAt(i, m)
      v.set(x, y + h, z)
      if (hidden) s.set(0.0001, 0.0001, 0.0001); else s.set(1, 1, 1)
      m.compose(v, q, s)
      this.caps.setMatrixAt(i, m)
    }
    this.columns.instanceMatrix.needsUpdate = true
    this.caps.instanceMatrix.needsUpdate = true
    if (this.columns.instanceColor) this.columns.instanceColor.needsUpdate = true
    if (this.caps.instanceColor) this.caps.instanceColor.needsUpdate = true
  }

  setBaseFromRelief(relief, exaggeration) {
    for (let i = 0; i < this.count; i++) {
      const s = this.sites[i]
      const h = relief ? relief.heightAt(s.x, s.y) : 0
      this.baseY[i] = h == null ? 0 : h * exaggeration * SCALE
    }
  }

  /**
   * Assign per-site values.
   * @param {Float32Array|null} values normalised 0..1 (null = uniform height)
   * @param {(i:number, v:number)=>THREE.Color} colorFn
   * @param {number} heightScale max column height in scene units
   * @param {Uint8Array|null} hiddenMask 1 = hide the site
   */
  apply(values, colorFn, heightScale, hiddenMask = null) {
    for (let i = 0; i < this.count; i++) {
      this.hidden[i] = hiddenMask ? hiddenMask[i] : 0
      const v = values ? values[i] : null
      this.heights[i] = v == null || Number.isNaN(v) ? 0.6 : 1 + v * heightScale
      const c = colorFn(i, v)
      this.setColor(i, c.r, c.g, c.b)
    }
    this._writeMatrices()
  }

  setVisible(v) {
    this.columns.visible = v
    this.caps.visible = v
  }

  setHighlight(index) {
    if (this.highlightIndex === index && index < 0) return
    const prev = this.highlightIndex
    this.highlightIndex = index
    const restore = (i) => {
      if (i < 0) return
      this._c.setRGB(this.colors[i * 3], this.colors[i * 3 + 1], this.colors[i * 3 + 2])
      this.columns.setColorAt(i, this._c)
      this.caps.setColorAt(i, this._c)
    }
    restore(prev)
    if (index >= 0) {
      this.columns.setColorAt(index, this.highlightColor)
      this.caps.setColorAt(index, this.highlightColor)
    }
    this.columns.instanceColor.needsUpdate = true
    this.caps.instanceColor.needsUpdate = true
  }

  siteIndexAt(intersection) {
    return intersection.instanceId ?? -1
  }

  dispose() {
    this.columns.geometry.dispose()
    this.columns.material.dispose()
    this.caps.geometry.dispose()
    this.caps.material.dispose()
  }
}
