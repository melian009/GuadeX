import { rampGradientCss } from '../lib/colors.js'

const $ = (sel, root = document) => root.querySelector(sel)

/**
 * Binds the static HTML skeleton to app callbacks and provides small
 * imperative helpers for updating dynamic pieces.
 */
export class AppUI {
  constructor(handlers = {}) {
    this.handlers = handlers
    this.el = {
      canvas: $('#scene'),
      loading: $('#loading'),
      loadingText: $('#loading-text'),
      panel: $('#panel'),
      layerToggles: $('#layer-toggles'),
      chkRelief: $('#chk-relief'),
      rangeExagg: $('#range-exagg'),
      outExagg: $('#out-exagg'),
      rangeHeight: $('#range-height'),
      outHeight: $('#out-height'),
      chkOnlyData: $('#chk-onlydata'),
      chkLabels: $('#chk-labels'),
      selSource: $('#sel-source'),
      selMetric: $('#sel-metric'),
      selRamp: $('#sel-ramp'),
      timeWrap: $('#time-wrap'),
      rangeTime: $('#range-time'),
      outTime: $('#out-time'),
      dropZone: $('#drop-zone'),
      fileResults: $('#file-results'),
      btnLoadFile: $('#btn-load-file'),
      btnDemoTime: $('#btn-load-demo-timeseries'),
      btnDemoMetrics: $('#btn-load-demo-metrics'),
      btnClearResults: $('#btn-clear-results'),
      resultsStatus: $('#results-status'),
      search: $('#search'),
      siteCodes: $('#site-codes'),
      legend: $('#legend'),
      legendTitle: $('#legend-title'),
      legendBar: $('#legend-bar'),
      legendMin: $('#legend-min'),
      legendMax: $('#legend-max'),
      legendCaption: $('#legend-caption'),
      tooltip: $('#tooltip'),
      info: $('#info'),
      infoTitle: $('#info-title'),
      infoSub: $('#info-sub'),
      infoBody: $('#info-body'),
      infoClose: $('#info-close'),
      stats: $('#stats'),
      btnReset: $('#btn-reset'),
      btnPanel: $('#btn-panel'),
    }
    this._bind()
  }

  _bind() {
    const h = this.handlers
    this.el.chkRelief.addEventListener('change', () => h.onReliefToggle?.(this.el.chkRelief.checked))
    this.el.rangeExagg.addEventListener('input', () => {
      const v = Number(this.el.rangeExagg.value)
      this.el.outExagg.textContent = `${v}×`
      h.onExaggeration?.(v)
    })
    this.el.rangeHeight.addEventListener('input', () => {
      const v = Number(this.el.rangeHeight.value)
      this.el.outHeight.textContent = String(v)
      h.onHeightScale?.(v)
    })
    this.el.chkOnlyData.addEventListener('change', () => h.onOnlyData?.(this.el.chkOnlyData.checked))
    this.el.chkLabels.addEventListener('change', () => h.onLabels?.(this.el.chkLabels.checked))
    this.el.selSource.addEventListener('change', () => h.onSourceChange?.(this.el.selSource.value))
    this.el.selMetric.addEventListener('change', () => h.onMetricChange?.(this.el.selMetric.value))
    this.el.selRamp.addEventListener('change', () => h.onRampChange?.(this.el.selRamp.value))
    this.el.rangeTime.addEventListener('input', () => h.onTimeChange?.(Number(this.el.rangeTime.value)))
    this.el.infoClose.addEventListener('click', () => this.hideInfo())
    this.el.btnReset.addEventListener('click', () => h.onResetView?.())
    this.el.btnPanel.addEventListener('click', () => this.togglePanel())
    this.el.btnClearResults.addEventListener('click', () => h.onClearResults?.())

    for (const tab of this.el.panel.querySelectorAll('.tab')) {
      tab.addEventListener('click', () => this.setTab(tab.dataset.tab))
    }

    this.el.btnLoadFile.addEventListener('click', () => this.el.fileResults.click())
    this.el.fileResults.addEventListener('change', (e) => {
      const file = e.target.files?.[0]
      if (file) h.onLoadFile?.(file)
      e.target.value = ''
    })
    this.el.btnDemoTime.addEventListener('click', () => h.onLoadDemo?.('timeseries'))
    this.el.btnDemoMetrics.addEventListener('click', () => h.onLoadDemo?.('metrics'))

    const dz = this.el.dropZone
    ;['dragenter', 'dragover'].forEach((ev) => dz.addEventListener(ev, (e) => { e.preventDefault(); dz.classList.add('drag') }))
    ;['dragleave', 'drop'].forEach((ev) => dz.addEventListener(ev, (e) => { e.preventDefault(); dz.classList.remove('drag') }))
    dz.addEventListener('drop', (e) => {
      const file = e.dataTransfer?.files?.[0]
      if (file) h.onLoadFile?.(file)
    })
    document.addEventListener('dragover', (e) => e.preventDefault())
    document.addEventListener('drop', (e) => e.preventDefault())

    this.el.search.addEventListener('change', () => h.onSearch?.(this.el.search.value))
    this.el.search.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') h.onSearch?.(this.el.search.value)
    })

    window.addEventListener('keydown', (e) => {
      if (e.target instanceof HTMLInputElement || e.target instanceof HTMLSelectElement) return
      if (e.key === 'r' || e.key === 'R') h.onResetView?.()
      if (e.key === 'p' || e.key === 'P') this.togglePanel()
      if (e.key === 'f' || e.key === 'F' || e.key === '/') { e.preventDefault(); this.el.search.focus() }
      if (e.key === 'Escape') this.hideInfo()
    })
  }

  setTab(name) {
    for (const tab of this.el.panel.querySelectorAll('.tab')) tab.classList.toggle('active', tab.dataset.tab === name)
    for (const panel of this.el.panel.querySelectorAll('.tab-panel')) panel.classList.toggle('active', panel.dataset.panel === name)
  }

  togglePanel() {
    this.el.panel.classList.toggle('collapsed')
  }

  setLoadingText(text) {
    this.el.loadingText.textContent = text
  }

  fadeLoading() {
    this.el.loading.classList.add('fade')
    setTimeout(() => this.el.loading.classList.add('hidden'), 600)
  }

  setLayerToggles(layers) {
    this.el.layerToggles.innerHTML = ''
    for (const layer of layers) {
      const label = document.createElement('label')
      label.className = 'toggle'
      label.innerHTML = `<input type="checkbox" ${layer.visible ? 'checked' : ''} />
        <span class="swatch" style="background:${layer.swatch ?? '#5aa9ff'}"></span>
        <span class="name">${layer.label}</span>
        <span class="count">${layer.count ?? ''}</span>`
      label.querySelector('input').addEventListener('change', (e) => this.handlers.onLayerToggle?.(layer.id, e.target.checked))
      this.el.layerToggles.appendChild(label)
    }
  }

  setMetricOptions(groups, selected) {
    const sel = this.el.selMetric
    sel.innerHTML = ''
    for (const group of groups) {
      if (!group.options.length) continue
      const og = document.createElement('optgroup')
      og.label = group.label
      for (const opt of group.options) {
        const o = document.createElement('option')
        o.value = opt.key
        o.textContent = opt.label
        if (opt.key === selected) o.selected = true
        og.appendChild(o)
      }
      sel.appendChild(og)
    }
  }

  setTime(steps, index) {
    if (!steps || steps.length <= 1) {
      this.el.timeWrap.classList.add('hidden')
      return
    }
    this.el.timeWrap.classList.remove('hidden')
    this.el.rangeTime.min = '0'
    this.el.rangeTime.max = String(steps.length - 1)
    this.el.rangeTime.value = String(index)
    this.el.outTime.textContent = steps[index] ?? '—'
  }

  setResultsStatus(text, kind = '') {
    this.el.resultsStatus.textContent = text
    this.el.resultsStatus.className = `status ${kind}`
  }

  setLegend({ title, ramp, min, max, caption, categorical, categories }) {
    if (!title) { this.el.legend.classList.add('hidden'); return }
    this.el.legend.classList.remove('hidden')
    this.el.legendTitle.textContent = title
    if (categorical && categories) {
      this.el.legendBar.style.background = `linear-gradient(90deg, ${categories.map((c) => c.color).join(', ')})`
      this.el.legendMin.textContent = ''
      this.el.legendMax.textContent = ''
      this.el.legendCaption.textContent = categories.map((c) => c.label).join(' · ')
    } else {
      this.el.legendBar.style.background = rampGradientCss(ramp)
      this.el.legendMin.textContent = min
      this.el.legendMax.textContent = max
      this.el.legendCaption.textContent = caption ?? ''
    }
  }

  setStats(items) {
    this.el.stats.innerHTML = items.map((i) => `<span>${i.label} <b>${i.value}</b></span>`).join('')
  }

  setSiteCodes(ids) {
    this.el.siteCodes.innerHTML = ids.slice(0, 4000).map((id) => `<option value="${id}"></option>`).join('')
  }

  setSearchValue(value) {
    this.el.search.value = value
  }

  showInfo({ title, subtitle, html }) {
    this.el.infoTitle.textContent = title
    this.el.infoSub.textContent = subtitle
    this.el.infoBody.innerHTML = html
    this.el.info.classList.remove('hidden')
  }

  hideInfo() {
    this.el.info.classList.add('hidden')
    this.handlers.onInfoClose?.()
  }

  showTooltip(html, clientX, clientY) {
    const t = this.el.tooltip
    t.innerHTML = html
    t.classList.remove('hidden')
    const pad = 14
    const rect = t.getBoundingClientRect()
    let x = clientX + pad
    let y = clientY + pad
    if (x + rect.width > window.innerWidth - 8) x = clientX - rect.width - pad
    if (y + rect.height > window.innerHeight - 8) y = clientY - rect.height - pad
    t.style.left = `${x}px`
    t.style.top = `${y}px`
  }

  hideTooltip() {
    this.el.tooltip.classList.add('hidden')
  }
}
