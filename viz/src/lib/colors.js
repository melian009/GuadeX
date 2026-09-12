import * as THREE from 'three'

/** Named colour ramps as lists of hex stops (low -> high). */
export const RAMPS = {
  viridis: ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'],
  turbo: ['#30123b', '#4145ab', '#4675ed', '#39a2fc', '#1bcfd4', '#24eca6', '#61fc6c', '#a4fc3b', '#d1e834', '#f3c63a', '#fe9b2d', '#f36315', '#d93806', '#b11901', '#7a0402'],
  magma: ['#000004', '#1c1044', '#4f127b', '#812581', '#b5367a', '#e55063', '#fb8761', '#fec287', '#fbfdbf'],
  rdbu: ['#313695', '#4575b4', '#74add1', '#abd9e9', '#e0f3f8', '#ffffbf', '#fee090', '#fdae61', '#f46d43', '#d73027', '#a50026'],
  terrain: ['#2e6f3e', '#5d9d54', '#9fc46b', '#d8d98a', '#c2a05f', '#9a7046', '#f2f2f2'],
  categorical: ['#4cc9f0', '#f72585', '#ffb703', '#2ec4b6', '#8e7dbe', '#ff6b6b', '#f9c74f', '#43aa8b', '#577590', '#f9844a', '#90be6d', '#c77dff'],
  native: ['#0b3d2e', '#14634b', '#2e9e74', '#6fd6ac', '#c9f5e4'],
  exotic: ['#3b1d00', '#7a3b00', '#c46a00', '#f0a63a', '#ffe1a8'],
}

const cache = new Map()

function hexToRgb(hex) {
  const n = parseInt(hex.slice(1), 16)
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255]
}

/** Interpolate a ramp, t in [0, 1]; returns [r, g, b] 0..255. */
export function rampRgb(name, t) {
  const stops = RAMPS[name] ?? RAMPS.viridis
  let table = cache.get(name)
  if (!table) {
    table = stops.map(hexToRgb)
    cache.set(name, table)
  }
  const x = Math.max(0, Math.min(1, t)) * (table.length - 1)
  const i = Math.floor(x)
  const f = x - i
  const a = table[i]
  const b = table[Math.min(i + 1, table.length - 1)]
  return [
    Math.round(a[0] + (b[0] - a[0]) * f),
    Math.round(a[1] + (b[1] - a[1]) * f),
    Math.round(a[2] + (b[2] - a[2]) * f),
  ]
}

export function rampCss(name, t) {
  const [r, g, b] = rampRgb(name, t)
  return `rgb(${r}, ${g}, ${b})`
}

/** CSS linear-gradient string for a ramp, used by the legend. */
export function rampGradientCss(name, steps = 12) {
  const parts = []
  for (let i = 0; i < steps; i++) {
    parts.push(`${rampCss(name, i / (steps - 1))} ${(i / (steps - 1) * 100).toFixed(1)}%`)
  }
  return `linear-gradient(90deg, ${parts.join(', ')})`
}

export function rampThreeColor(name, t, target = new THREE.Color()) {
  const [r, g, b] = rampRgb(name, t)
  return target.setRGB(r / 255, g / 255, b / 255)
}

/** Deterministic categorical colour from an arbitrary key. */
export function categoricalColor(key) {
  const stops = RAMPS.categorical
  let h = 2166136261
  const s = String(key ?? '')
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i)
    h = Math.imul(h, 16777619)
  }
  return stops[Math.abs(h) % stops.length]
}
