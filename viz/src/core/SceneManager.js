import * as THREE from 'three'
import { OrbitControls } from 'three/addons/controls/OrbitControls.js'
import { sceneCentre, sceneExtent } from '../lib/geo.js'

const BG = 0x070b12

/**
 * Owns the renderer, camera, controls, lighting and picking.
 * Layers register themselves as pickable objects and receive hover/click callbacks.
 */
export class SceneManager {
  constructor(canvas) {
    this.canvas = canvas
    this.pickables = new Map()
    this.frameCallbacks = new Set()

    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, powerPreference: 'high-performance' })
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2))
    this.renderer.outputColorSpace = THREE.SRGBColorSpace
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping
    this.renderer.toneMappingExposure = 1.08

    this.scene = new THREE.Scene()
    this.scene.background = new THREE.Color(BG)
    this.scene.fog = new THREE.FogExp2(BG, 0.00055)

    this.camera = new THREE.PerspectiveCamera(46, 1, 0.5, 8000)
    this.camera.position.set(0, 320, 480)

    this.controls = new OrbitControls(this.camera, canvas)
    this.controls.enableDamping = true
    this.controls.dampingFactor = 0.075
    this.controls.rotateSpeed = 0.55
    this.controls.zoomSpeed = 0.9
    this.controls.panSpeed = 0.7
    this.controls.maxPolarAngle = Math.PI * 0.495
    this.controls.minDistance = 6
    this.controls.maxDistance = 2400

    this._addEnvironment()
    this._bindPointer()
    this._bindResize()

    this._clock = new THREE.Clock()
    this._loop = this._loop.bind(this)
    this.renderer.setAnimationLoop(this._loop)
  }

  _addEnvironment() {
    const hemi = new THREE.HemisphereLight(0xcfe0ff, 0x141c2b, 0.85)
    this.scene.add(hemi)

    const sun = new THREE.DirectionalLight(0xffffff, 1.65)
    sun.position.set(0.7, 1.1, 0.45)
    this.scene.add(sun)

    const fill = new THREE.DirectionalLight(0x7fb2ff, 0.35)
    fill.position.set(-0.8, 0.5, -0.6)
    this.scene.add(fill)

    this.groundMaterial = new THREE.MeshStandardMaterial({
      color: 0x0b1220,
      roughness: 1,
      metalness: 0,
      transparent: true,
      opacity: 0.96,
    })
    this.ground = new THREE.Mesh(new THREE.PlaneGeometry(6000, 6000), this.groundMaterial)
    this.ground.rotation.x = -Math.PI / 2
    this.ground.position.y = -0.35
    this.scene.add(this.ground)

    this.grid = new THREE.GridHelper(2000, 40, 0x1d2c44, 0x142033)
    this.grid.position.y = -0.3
    this.grid.material.transparent = true
    this.grid.material.opacity = 0.5
    this.scene.add(this.grid)
  }

  _bindResize() {
    const onResize = () => {
      const w = this.canvas.clientWidth || window.innerWidth
      const h = this.canvas.clientHeight || window.innerHeight
      this.renderer.setSize(w, h, false)
      this.camera.aspect = w / h
      this.camera.updateProjectionMatrix()
      this.dispatchResize(w, h)
    }
    this._onResize = onResize
    window.addEventListener('resize', onResize)
    if (typeof ResizeObserver !== 'undefined') {
      this._ro = new ResizeObserver(onResize)
      this._ro.observe(this.canvas)
    }
    onResize()
  }

  onResize(cb) {
    this.resizeCallbacks = this.resizeCallbacks || new Set()
    this.resizeCallbacks.add(cb)
    cb(this.canvas.clientWidth, this.canvas.clientHeight)
  }

  dispatchResize(w, h) {
    if (this.resizeCallbacks) for (const cb of this.resizeCallbacks) cb(w, h)
  }

  _bindPointer() {
    const ndc = new THREE.Vector2()
    const raycaster = new THREE.Raycaster()
    this._raycaster = raycaster

    const setNdc = (event) => {
      const rect = this.canvas.getBoundingClientRect()
      ndc.x = ((event.clientX - rect.left) / rect.width) * 2 - 1
      ndc.y = -((event.clientY - rect.top) / rect.height) * 2 + 1
    }

    const hitTest = () => {
      raycaster.setFromCamera(ndc, this.camera)
      const objects = [...this.pickables.keys()].filter((o) => o.visible)
      const hits = raycaster.intersectObjects(objects, false)
      return hits.find((h) => this.pickables.get(h.object))
    }

    let downX = 0, downY = 0, moved = false
    this.canvas.addEventListener('pointerdown', (e) => {
      downX = e.clientX; downY = e.clientY; moved = false
    })
    this.canvas.addEventListener('pointermove', (e) => {
      if (Math.abs(e.clientX - downX) > 4 || Math.abs(e.clientY - downY) > 4) moved = true
      setNdc(e)
      const hit = hitTest()
      const entry = hit ? this.pickables.get(hit.object) : null
      const key = hit ? (entry?.hitKey?.(hit) ?? `${hit.object.uuid}:${hit.instanceId ?? ''}`) : null
      if (this._hoverKey !== key) {
        this._hoverKey = key
        this._currentHoverEntry?.onHover?.(null, null)
        this._currentHoverEntry = entry
        entry?.onHover?.(hit, e)
        this.canvas.style.cursor = hit ? 'pointer' : 'default'
      }
      if (entry?.onMove) entry.onMove(hit, e)
    })
    this.canvas.addEventListener('pointerleave', () => {
      this._currentHoverEntry?.onHover?.(null, null)
      this._currentHoverEntry = null
      this._hoverKey = null
      this.canvas.style.cursor = 'default'
    })
    this.canvas.addEventListener('pointerup', (e) => {
      if (moved) return
      setNdc(e)
      const hit = hitTest()
      if (hit) this.pickables.get(hit.object)?.onClick?.(hit, e)
    })
  }

  registerPickable(object, handlers) {
    this.pickables.set(object, handlers)
  }

  unregisterPickable(object) {
    this.pickables.delete(object)
  }

  add(object) {
    this.scene.add(object)
  }

  remove(object) {
    this.scene.remove(object)
  }

  onFrame(cb) {
    this.frameCallbacks.add(cb)
  }

  fitBounds(bounds, { padding = 1.25 } = {}) {
    const [cx, cz] = sceneCentre(bounds)
    const { width, depth } = sceneExtent(bounds)
    const radius = Math.max(width, depth) / 2
    const fov = THREE.MathUtils.degToRad(this.camera.fov)
    const dist = (radius / Math.tan(fov / 2)) * padding
    this.controls.target.set(cx, 0, cz)
    this.camera.position.set(cx, dist * 0.85, cz + dist * 0.95)
    this.controls.update()
    this._home = {
      pos: this.camera.position.clone(),
      target: this.controls.target.clone(),
    }
  }

  focus(point, distance = 70) {
    const target = new THREE.Vector3(point.x, point.y, point.z)
    const dir = this.camera.position.clone().sub(this.controls.target).normalize()
    this.controls.target.copy(target)
    this.camera.position.copy(target.clone().add(dir.multiplyScalar(distance)))
    this.controls.update()
  }

  resetView() {
    if (!this._home) return
    this.camera.position.copy(this._home.pos)
    this.controls.target.copy(this._home.target)
    this.controls.update()
  }

  projectToScreen(vec3) {
    const v = vec3.clone().project(this.camera)
    const rect = this.canvas.getBoundingClientRect()
    return {
      x: (v.x * 0.5 + 0.5) * rect.width,
      y: (-v.y * 0.5 + 0.5) * rect.height,
      visible: v.z < 1,
    }
  }

  _loop() {
    const dt = this._clock.getDelta()
    this.controls.update()
    for (const cb of this.frameCallbacks) cb(dt)
    this.renderer.render(this.scene, this.camera)
  }

  dispose() {
    this.renderer.setAnimationLoop(null)
    window.removeEventListener('resize', this._onResize)
    this._ro?.disconnect()
    this.controls.dispose()
    this.renderer.dispose()
  }
}
