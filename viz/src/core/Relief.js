import * as THREE from 'three'
import { sceneX, sceneZ, SCALE } from '../lib/geo.js'
import { rampThreeColor } from '../lib/colors.js'

const SENTINEL = -32768

function decodeInt16(base64) {
  const bin = atob(base64)
  const bytes = new Uint8Array(bin.length)
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i)
  return new Int16Array(bytes.buffer)
}

/**
 * Interpolated relief surface sampled on a regular grid and masked to the
 * catchment union. Values are metres; -32768 marks cells outside the basin.
 */
export class Relief {
  constructor(json) {
    this.width = json.width
    this.height = json.height
    this.cellSize = json.cellSize
    this.minX = json.origin[0]
    this.minY = json.origin[1]
    this.data = decodeInt16(json.data)

    let min = Infinity, max = -Infinity
    for (const v of this.data) {
      if (v === SENTINEL) continue
      if (v < min) min = v
      if (v > max) max = v
    }
    this.minElev = min
    this.maxElev = max
    this.exaggeration = 1
  }

  isMasked(i, j) {
    if (i < 0 || j < 0 || i >= this.width || j >= this.height) return true
    return this.data[j * this.width + i] === SENTINEL
  }

  /** Bilinear sampled elevation in metres, or null outside the basin. */
  heightAt(xMetres, yMetres) {
    const fi = (xMetres - this.minX) / this.cellSize
    const fj = (yMetres - this.minY) / this.cellSize
    const i0 = Math.floor(fi)
    const j0 = Math.floor(fj)
    const d = this.data
    const w = this.width
    const at = (i, j) => (i < 0 || j < 0 || i >= w || j >= this.height ? SENTINEL : d[j * w + i])
    const v00 = at(i0, j0)
    const v10 = at(i0 + 1, j0)
    const v01 = at(i0, j0 + 1)
    const v11 = at(i0 + 1, j0 + 1)
    const vals = [v00, v10, v01, v11].filter((v) => v !== SENTINEL)
    if (vals.length === 0) return null
    if (v00 !== SENTINEL && v10 !== SENTINEL && v01 !== SENTINEL && v11 !== SENTINEL) {
      const tx = fi - i0
      const ty = fj - j0
      const a = v00 + (v10 - v00) * tx
      const b = v01 + (v11 - v01) * tx
      return a + (b - a) * ty
    }
    return vals.reduce((s, v) => s + v, 0) / vals.length
  }

  buildMesh() {
    const { width, height, cellSize } = this
    const n = width * height
    const positions = new Float32Array(n * 3)
    const baseElev = new Float32Array(n)
    const colors = new Float32Array(n * 3)
    const range = Math.max(1, this.maxElev - this.minElev)
    const color = new THREE.Color()

    for (let j = 0; j < height; j++) {
      const y = this.minY + j * cellSize
      for (let i = 0; i < width; i++) {
        const idx = j * width + i
        const e = this.data[idx]
        const p = idx * 3
        positions[p] = sceneX(this.minX + i * cellSize)
        positions[p + 1] = 0
        positions[p + 2] = sceneZ(y)
        baseElev[idx] = e === SENTINEL ? -9999 : e
        if (e !== SENTINEL) {
          rampThreeColor('terrain', (e - this.minElev) / range, color)
          colors[p] = color.r; colors[p + 1] = color.g; colors[p + 2] = color.b
        }
      }
    }

    const indices = []
    for (let j = 0; j < height - 1; j++) {
      for (let i = 0; i < width - 1; i++) {
        const a = j * width + i
        const b = a + 1
        const c = a + width
        const d = c + 1
        if (this.isMasked(i, j) || this.isMasked(i + 1, j) || this.isMasked(i, j + 1) || this.isMasked(i + 1, j + 1)) continue
        indices.push(a, c, b, b, c, d)
      }
    }

    const geometry = new THREE.BufferGeometry()
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3))
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3))
    geometry.setIndex(indices)

    const material = new THREE.MeshStandardMaterial({
      vertexColors: true,
      roughness: 0.98,
      metalness: 0,
      flatShading: false,
      side: THREE.DoubleSide,
    })

    this.mesh = new THREE.Mesh(geometry, material)
    this.mesh.name = 'relief'
    this.baseElev = baseElev
    this.applyExaggeration(this.exaggeration)
    return this.mesh
  }

  applyExaggeration(exaggeration) {
    this.exaggeration = exaggeration
    if (!this.mesh) return
    const pos = this.mesh.geometry.attributes.position
    const arr = pos.array
    const k = exaggeration * SCALE
    for (let i = 0; i < this.baseElev.length; i++) {
      const e = this.baseElev[i]
      arr[i * 3 + 1] = e <= -9998 ? 0 : e * k
    }
    pos.needsUpdate = true
    this.mesh.geometry.computeVertexNormals()
    this.mesh.geometry.computeBoundingSphere()
  }

  setVisible(v) {
    if (this.mesh) this.mesh.visible = v
  }

  /** Scene-space Y for a geographic point (scene units). */
  sceneHeightAt(xMetres, yMetres, exaggeration = this.exaggeration) {
    const h = this.heightAt(xMetres, yMetres)
    if (h == null) return 0
    return h * exaggeration * SCALE
  }
}
