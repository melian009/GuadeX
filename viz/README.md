# GuadeX · Guadalquivir Basin 3-D Site Explorer

An interactive, browser-based 3-D view of the whole Guadalquivir river basin and the
**1,037 GuadeX sampling sites**, built so that *any* per-site data — including future
simulation outputs — can be attached, explored and compared by users.

This lives in its own directory (`viz/`) and is independent of the Julia analysis code.

---

## Quick start

```bash
cd viz
npm install
npm run data      # converts the shapefiles in ../data into ./public/data (already run)
npm run dev       # opens a dev server on http://localhost:5173
```

Production build / preview:

```bash
npm run build
npm run preview
```

---

## Hosting on GitHub Pages (no server needed)

The app is fully static, so GitHub Pages can host it for free. A workflow is provided
at `.github/workflows/deploy-viz.yml`; it installs dependencies, regenerates the web data
from the tracked GIS sources, builds the site, and deploys `viz/dist`.

One-time setup (repo **Settings**):

1. **Pages → Build and deployment → Source = “GitHub Actions”**.
2. **Actions → General → Workflow permissions**: the workflow declares its own
   `pages: write` / `id-token: write` permissions, so the default read-only setting is fine.
   Just make sure Actions are enabled.
3. The repository must be **public** (Pages on a private repo requires a paid plan).
4. Push to `master` (or run the workflow manually). The site appears at
   `https://<owner>.github.io/GuadeX/`.

`vite.config.js` uses `base: './'` and the loader uses relative paths, so the project-page
subpath works without further configuration. Because the workflow runs `npm run data`, the
generated `public/data/` does not have to be committed (you may add it to `.gitignore`).

### Can users upload their own simulation results?

Yes — entirely in the browser, with no server involved:

* The **Data → Simulation results** panel accepts a JSON/CSV file (choose or drag-and-drop).
* Uploaded files never leave the user's machine; the static host only serves the app.
* Uploads are per-session: refreshing clears them. To share a result set, commit the JSON
  next to the app (e.g. `public/data/my_results.json`) and send a deep link such as
  `?results=./data/my_results.json&site=1.1.2`, or host the JSON anywhere with CORS enabled
  and pass its URL to `?results=` / `window.GuadeX.loadResults(url)`.

Requirements for a user file: it must be keyed by the site code (`CODIGO`, e.g. `"1.1.2"`);
codes that do not match one of the 1,037 sites are ignored. Persisting uploads for other
users centrally would require a backend, which Pages does not provide.

---

## What is on screen

* **Catchments** (`Cuencas_masas_agua_4c`, 468 polygons) — drape the basin and act as a
  choropleth of the mean value of the active variable per water body.
* **Sampling sites** (`puntos muestreo con id masa`, 1,037 points) — instanced 3-D columns.
  Height and colour encode a selectable variable; clicking opens the full attribute card.
* **Rivers / streams** (`SW_Line_4C_`, 360 lines), **reservoirs & lakes** (`SWB_Lago`, 97),
  **groundwater bodies** (`GWB_4C`, 90), **transitional** (13) and **coastal** (3 ES050) waters.
* **Interpolated relief** — an IDW surface built from the 776 sites that carry an altitude,
  masked to the catchment union. It is a reading aid, **not a survey DEM**.

Every site carries the full attribute set from the GIS archive (species presence, species
densities, richness, ZIC/BCC/VA indices, habitat and water-quality fields), shown in the
site panel and available as colour/height variables.

---

## Connecting simulation results

Results are loaded at runtime either from a file, from a URL, or programmatically.
Three shapes are supported and normalised into one model.

**1 · Time series (JSON)** — most useful for scenario/time outputs:

```json
{
  "name": "Native extinction risk",
  "unit": "probability",
  "steps": ["2026", "2050", "2100"],
  "data": {
    "1.1.2": [0.11, 0.24, 0.61],
    "1.1.3": [0.08, 0.19, 0.55]
  }
}
```

A time slider appears; picking a step re-colours and re-scales the sites.

**2 · Per-site metrics (JSON)** — any key/value object per site:

```json
{
  "name": "Simulation output",
  "data": {
    "1.1.2": { "risk": 0.31, "population": 4200, "status": "at-risk" },
    "1.1.3": { "risk": 0.22, "population": 9800, "status": "stable" }
  }
}
```

Each key becomes a selectable variable. Non-numeric keys are treated as categories.

**3 · CSV** — long (`CODIGO,step,value`) or wide (`CODIGO,metricA,metricB`) format.
The site-id column may be named `CODIGO`, `id`, `site`, `code`, …

Values may also be scalar numbers; a single `value` key is created.

### Programmatic API

```js
window.GuadeX.setResults(jsonObjectOrCsvString)
await window.GuadeX.loadResults('./data/my_results.json')
window.GuadeX.selectSite('1.1.2')
window.GuadeX.setMetric('extinction_risk', 'results')
window.GuadeX.setRamp('turbo')
window.GuadeX.clearResults()
window.GuadeX.store.sites        // full site array
window.GuadeX.viewer             // SceneManager (Three.js)
```

### Shareable deep links

```
?results=./data/results.demo-timeseries.json&metric=extinction_risk&site=1.1.2&ramp=turbo
```

`GUADEX_READY` is announced via `window.dispatchEvent(new CustomEvent('guadex:ready'))`.

Two **clearly labelled synthetic demo files** ship so the pipeline can be tested
(`public/data/results.demo-timeseries.json`, `results.demo-metrics.json`). They are
generated from exotic richness, elevation and a hash of the site code — replace them.

---

## Data pipeline

`scripts/build-data.mjs` reads the raw shapefiles and tables and writes compact JSON to
`public/data/`:

| Source (relative to `../data/`) | Output | Notes |
| :--- | :--- | :--- |
| `GIS/capas GIS/capa_puntos muestreo con id masa/…shp` | `sites.json` | 1,037 sites + all attributes |
| `Version_02-01-2026-Masas/Cuencas_masas_agua_4c.shp` | `catchments.json` | 468 catchment polygons |
| `Version_02-01-2026-Masas/GWB_4C.shp` | `gwb.json` | 90 groundwater bodies |
| `Version_02-01-2026-Masas/SW_Line_4C_.shp` | `rivers.json` | 360 river lines |
| `Version_02-01-2026-Masas/SWB_Lago.shp` | `reservoirs.json` | 97 lakes/reservoirs |
| `Version_02-01-2026-Masas/SWB_Transicion.shp` | `transitional.json` | 13 transitional waters |
| `Version_02-01-2026-Masas/SWB_Costera.shp` | `coastal.json` | ES050 coastal waters |
| `ConnectivityUTM.csv` | (joined into `sites.json`) | supplies site altitude |
| `GIS/info atributos capas GIS/info todas-449 masas aguas.xlsx` | `waterbodies.json` | water-body names/areas |

Details:

* CRS is **EPSG:25830 (ETRS89 / UTM zone 30N, metres)**. Coordinates are recentred around the
  basin centre (`manifest.origin`) and rounded to whole metres so the browser never handles
  the large UTM magnitudes.
* Geometry is simplified with Douglas–Peucker (per-layer tolerances in `TOLERANCE`).
* Relief is generated by rasterising the catchment union into a 1 km mask, then IDW
  (k = 10, 1/d³) interpolating the 776 site altitudes, stored as a base64 `Int16` grid.
* Text fields are repaired from the shapefile's latin1/UTF-8 mis-decoding.

Re-run any time with `npm run data`. `npm run smoke` runs headless logic checks over the
generated relief and the results parser.

---

## Architecture

```
viz/
  index.html
  src/
    main.js                 app controller + state
    core/
      SceneManager.js       renderer, camera, lights, picking
      Relief.js             masked IDW heightfield surface
    layers/
      PolygonLayer.js       merged, per-feature-colourable, pickable polygons
      LineLayer.js          screen-width river lines
      SitesLayer.js         instanced site columns
    data/
      DataStore.js          dataset loading/lookup
      ResultsModel.js       results normalisation (JSON/CSV → uniform model)
    lib/
      geo.js                coordinate transforms
      colors.js             colour ramps
    ui/AppUI.js             DOM bindings + panels
  scripts/
    build-data.mjs          GIS → web data
    inspect-data.mjs        attribute-table inspector (debug)
    smoke-test.mjs          logic checks
  public/data/              generated data (safe to regenerate)
```

### Known caveats

* Site→catchment linkage matches `ID_masa`/`EUMASCod` for **999 of 1,037** sites; the rest
  sit in water bodies outside the ES050 catchment set and are still shown as points.
* The relief is interpolated from sparse site altitudes — it conveys the landscape but is
  not terrain data. Turn it off in **Display → Interpolated relief** to work on a flat plane.
* The shipped demo results are synthetic and labelled as such.
