import * as THREE from 'three'
import { sceneX, sceneZ } from '../lib/geo.js'

function toVec2s(ring) {
  const out = new Array(ring.length)
  for (let i = 0; i < ring.length; i++) out[i] = new THREE.Vector2(ring[i][0], ring[i][1])
  return out
}

function iteratePolygons(geometry) {
  if (!geometry) return []
  if (geometry.type === 'Polygon') return [geometry.coordinates]
  if (geometry.type === 'MultiPolygon') return geometry.coordinates
  return []
}

/**
 * A merged, pickable, per-feature-colourable polygon layer.
 * Coordinates are recentred metres; height is draped via `heightAt(x, y)`.
 */
export class PolygonLayer {
  constructor({
    name,
    features,
    color = 0x2a4a6b,
    opacity = 0.88,
    heightAt = null,
    outline = true,
    outlineColor = 0x6f9dcf,
    outlineOpacity = 0.4,
    outlineWidth = 1,
    renderOrder = 0,
  }) {
    this.name = name
    this.features = features
    this.heightAt = heightAt
    this.baseColor = new THREE.Color(color)
    this.opacity = opacity
    this.renderOrder = renderOrder

    this.featureIds = []
    this.featureProps = []
    this.vertexRanges = []
    this.faceFeature = []
    this.centroids = []
    this.featureColors = null
    this.colorFn = null
    this.highlightIndex = -1
    this._highlightColor = new THREE.Color(0xffffff)

    this._buildMesh()
    if (outline) this._buildOutline(outlineColor, outlineOpacity, outlineWidth)
    this._buildPickingMesh()
  }

  get mesh() { return this.mesh_ }
  get outlineMesh() { return this.outline_ }
  get pickTarget() { return this.pickMesh_ }

  _buildMesh() {
    const positions = []
    const geo = []
    const colors = []
    const indices = []
    const faceFeature = []

    for (let fi = 0; fi < this.features.length; fi++) {
      const feature = this.features[fi]
      const id = feature.properties?.id ?? String(fi)
      this.featureIds.push(id)
      this.featureProps.push(feature.properties ?? {})
      const vertexStart = positions.length / 3
      let centroidX = 0, centroidY = 0, centroidN = 0

      for (const rings of iteratePolygons(feature.geometry)) {
        if (!rings.length) continue
        const contour = toVec2s(rings[0])
        const holes = rings.slice(1).map(toVec2s)
        let tris
        try {
          tris = THREE.ShapeUtils.triangulateShape(contour, holes)
        } catch {
          tris = null
        }
        const pts = contour.concat(...holes)
        const offset = positions.length / 3
        for (const p of pts) {
          const x = p.x, y = p.y
          const h = this.heightAt ? this.heightAt(x, y) : 0
          positions.push(sceneX(x), h, sceneZ(y))
          geo.push(x, y)
          colors.push(this.baseColor.r, this.baseColor.g, this.baseColor.b)
          centroidX += x; centroidY += y; centroidN++
        }
        if (tris) {
          for (const t of tris) {
            indices.push(t[0] + offset, t[1] + offset, t[2] + offset)
            faceFeature.push(fi)
          }
        }
      }
      this.vertexRanges.push([vertexStart, positions.length / 3 - vertexStart])
      this.centroids.push(
        centroidN
          ? new THREE.Vector3(sceneX(centroidX / centroidN), 0, sceneZ(centroidY / centroidN))
          : new THREE.Vector3(),
      )
    }

    const geometry = new THREE.BufferGeometry()
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3))
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3))
    geometry.setIndex(indices)
    geometry.computeVertexNormals()
    geometry.computeBoundingSphere()

    const material = new THREE.MeshStandardMaterial({
      vertexColors: true,
      transparent: this.opacity < 1,
      opacity: this.opacity,
      roughness: 0.94,
      metalness: 0,
      side: THREE.DoubleSide,
      depthWrite: this.opacity >= 1,
      polygonOffset: true,
      polygonOffsetFactor: 1,
      polygonOffsetUnits: 1,
    })

    this.mesh_ = new THREE.Mesh(geometry, material)
    this.mesh_.name = this.name
    this.mesh_.renderOrder = this.renderOrder
    this.vertexGeo = new Float32Array(geo)
    this.faceFeature = faceFeature
    this.featureColors = new Float32Array(this.features.length * 3)
    for (let i = 0; i < this.features.length; i++) {
      this.featureColors[i * 3] = this.baseColor.r
      this.featureColors[i * 3 + 1] = this.baseColor.g
      this.featureColors[i * 3 + 2] = this.baseColor.b
    }
  }

  _buildOutline(color, opacity, width) {
    const positions = []
    const geo = []
    for (const feature of this.features) {
      for (const rings of iteratePolygons(feature.geometry)) {
        for (const ring of rings) {
          for (let i = 0; i < ring.length - 1; i++) {
            const a = ring[i]
            const b = ring[i + 1]
            const ha = (this.heightAt ? this.heightAt(a[0], a[1]) : 0) + 0.06
            const hb = (this.heightAt ? this.heightAt(b[0], b[1]) : 0) + 0.06
            positions.push(sceneX(a[0]), ha, sceneZ(a[1]))
            positions.push(sceneX(b[0]), hb, sceneZ(b[1]))
            geo.push(a[0], a[1], b[0], b[1])
          }
        }
      }
    }
    const geometry = new THREE.BufferGeometry()
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3))
    const material = new THREE.LineBasicMaterial({ color, transparent: opacity < 1, opacity, linewidth: width })
    this.outline_ = new THREE.LineSegments(geometry, material)
    this.outline_.name = `${this.name}-outline`
    this.outline_.renderOrder = this.renderOrder + 0.5
    this.outlineGeo = new Float32Array(geo)
  }

  /** Invisible triangles slightly above the fill, used for reliable picking. */
  _buildPickingMesh() {
    const geometry = this.mesh_.geometry.clone()
    const material = new THREE.MeshBasicMaterial({
      transparent: true,
      opacity: 0,
      depthWrite: false,
      side: THREE.DoubleSide,
    })
    this.pickMesh_ = new THREE.Mesh(geometry, material)
    this.pickMesh_.name = `${this.name}-pick`
    this.pickMesh_.userData.layer = this
  }

  /** Set a colour function: fn(featureIndex, properties) -> [r,g,b] in 0..1. */
  setColorFunction(fn) {
    this.colorFn = fn
    this.applyColors()
  }

  setBaseColor(color) {
    this.baseColor.set(color)
    if (!this.colorFn) {
      for (let i = 0; i < this.features.length; i++) {
        this.featureColors[i * 3] = this.baseColor.r
        this.featureColors[i * 3 + 1] = this.baseColor.g
        this.featureColors[i * 3 + 2] = this.baseColor.b
      }
      this.applyColors()
    }
  }

  setHighlight(index, color = 0xffffff) {
    this.highlightIndex = index
    this._highlightColor.set(color)
    this.applyColors()
  }

  applyColors() {
    if (!this.mesh_) return
    const attr = this.mesh_.geometry.attributes.color
    const arr = attr.array
    for (let fi = 0; fi < this.features.length; fi++) {
      let r, g, b
      if (fi === this.highlightIndex) {
        r = this._highlightColor.r; g = this._highlightColor.g; b = this._highlightColor.b
      } else if (this.colorFn) {
        const c = this.colorFn(fi, this.featureProps[fi])
        if (c.isColor) { r = c.r; g = c.g; b = c.b }
        else { r = c[0]; g = c[1]; b = c[2] }
      } else {
        r = this.featureColors[fi * 3]
        g = this.featureColors[fi * 3 + 1]
        b = this.featureColors[fi * 3 + 2]
      }
      const [start, count] = this.vertexRanges[fi]
      for (let i = start; i < start + count; i++) {
        arr[i * 3] = r; arr[i * 3 + 1] = g; arr[i * 3 + 2] = b
      }
    }
    attr.needsUpdate = true
  }

  /** Index of the feature hit by a raycast intersection, or -1. */
  featureIndexAt(intersection) {
    if (intersection.faceIndex == null) return -1
    return this.faceFeature[intersection.faceIndex] ?? -1
  }

  getFeatureId(index) {
    return this.featureIds[index]
  }

  setVisible(visible) {
    this.mesh_.visible = visible
    if (this.outline_) this.outline_.visible = visible
    this.pickMesh_.visible = visible
  }

  /** Re-drape every vertex after the vertical exaggeration changes. */
  refreshHeights() {
    if (!this.heightAt || !this.vertexGeo) return
    const meshes = [this.mesh_, this.pickMesh_]
    for (const mesh of meshes) {
      const arr = mesh.geometry.attributes.position.array
      for (let i = 0; i < this.vertexGeo.length / 2; i++) {
        arr[i * 3 + 1] = this.heightAt(this.vertexGeo[i * 2], this.vertexGeo[i * 2 + 1])
      }
      mesh.geometry.attributes.position.needsUpdate = true
      mesh.geometry.computeVertexNormals()
      mesh.geometry.computeBoundingSphere()
    }
    if (!this.outline_ || !this.outlineGeo) return
    const outArr = this.outline_.geometry.attributes.position.array
    for (let i = 0; i < this.outlineGeo.length / 2; i++) {
      outArr[i * 3 + 1] = (this.heightAt(this.outlineGeo[i * 2], this.outlineGeo[i * 2 + 1]) || 0) + 0.06
    }
    this.outline_.geometry.attributes.position.needsUpdate = true
    this.outline_.geometry.computeBoundingSphere()
  }

  dispose() {
    this.mesh_.geometry.dispose()
    this.mesh_.material.dispose()
    this.outline_?.geometry.dispose()
    this.outline_?.material.dispose()
    this.pickMesh_.geometry.dispose()
    this.pickMesh_.material.dispose()
  }
}
