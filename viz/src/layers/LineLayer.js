import * as THREE from 'three'
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js'
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js'
import { LineMaterial } from 'three/addons/lines/LineMaterial.js'
import { sceneX, sceneZ } from '../lib/geo.js'

function iterateLines(geometry) {
  if (!geometry) return []
  if (geometry.type === 'LineString') return [geometry.coordinates]
  if (geometry.type === 'MultiLineString') return geometry.coordinates
  return []
}

/** Screen-space-width line layer for the river network. */
export class LineLayer {
  constructor({
    name,
    features,
    color = 0x66c2ff,
    width = 1.6,
    opacity = 0.85,
    heightAt = null,
    renderOrder = 1,
  }) {
    this.name = name
    this.features = features
    this.heightAt = heightAt
    this.color = new THREE.Color(color)

    const positions = []
    this.segments = []
    for (const feature of features) {
      for (const line of iterateLines(feature.geometry)) {
        const geo = []
        for (let i = 0; i < line.length; i++) {
          const [x, y] = line[i]
          geo.push(x, y)
        }
        this.segments.push(geo)
        for (let i = 0; i < line.length - 1; i++) {
          const a = line[i]
          const b = line[i + 1]
          const ha = (heightAt ? heightAt(a[0], a[1]) : 0) + 0.14
          const hb = (heightAt ? heightAt(b[0], b[1]) : 0) + 0.14
          positions.push(sceneX(a[0]), ha, sceneZ(a[1]))
          positions.push(sceneX(b[0]), hb, sceneZ(b[1]))
        }
      }
    }

    const geometry = new LineSegmentsGeometry()
    geometry.setPositions(positions)

    this.material = new LineMaterial({
      color,
      linewidth: width,
      transparent: opacity < 1,
      opacity,
      worldUnits: false,
      dashed: false,
      alphaToCoverage: true,
    })
    this.line = new LineSegments2(geometry, this.material)
    this.line.name = name
    this.line.renderOrder = renderOrder
    this.line.computeLineDistances()
  }

  setResolution(width, height) {
    this.material.resolution.set(width, height)
  }

  setVisible(v) {
    this.line.visible = v
  }

  refreshHeights() {
    if (!this.heightAt || !this.segments) return
    const positions = []
    for (const line of this.segments) {
      for (let i = 0; i < line.length / 2 - 1; i++) {
        const ax = line[i * 2], ay = line[i * 2 + 1]
        const bx = line[i * 2 + 2], by = line[i * 2 + 3]
        const ha = (this.heightAt(ax, ay) || 0) + 0.14
        const hb = (this.heightAt(bx, by) || 0) + 0.14
        positions.push(sceneX(ax), ha, sceneZ(ay))
        positions.push(sceneX(bx), hb, sceneZ(by))
      }
    }
    this.line.geometry.setPositions(positions)
    this.line.computeLineDistances()
  }

  setColor(color) {
    this.material.color.set(color)
  }

  dispose() {
    this.line.geometry.dispose()
    this.material.dispose()
  }
}
