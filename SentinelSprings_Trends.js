/***********************************************************************
 * VEGETATION TREND ANALYSIS — Mann-Kendall & Sen's slope
 * Sentinel-2 SR (BOA), seasonal 90th-percentile composites, 2017-today
 *
 * Supported indices:
 *   NDVI  = (B8 - B4)  / (B8 + B4)    greenness
 *   NDMI  = (B8 - B11) / (B8 + B11)   vegetation water content
 *
 * When multiple indices are requested, MK and Sen's slope are computed
 * INDEPENDENTLY on each one (they are different physical quantities
 * with different statistical behaviour) and exported as SEPARATE
 * GeoTIFFs / assets, one per index.
 *
 * Workflow:
 *   1. User parameters             (the ONLY block each partner edits)
 *   2. AOI loading + display
 *   3. Sentinel-2 collection       + Cloud Score+ masking + index
 *   4. Annual seasonal p90 composites + per-pixel QC
 *   5. Time-series chart per AOI feature  (sanity check)
 *   6. Coverage diagnostic
 *   7. Mann-Kendall + Sen's slope
 *   8. Results & display
 *   9. Exports                     (results + annual composites)
 *
 * Notes:
 *   - Fully server-side / lazy.
 *   - Northern-hemisphere AOIs (season does not cross 1 Jan).
 ***********************************************************************/

/* ======================= 1. USER PARAMETERS ======================== */

// --- Area of Interest -------------------------------------------------
// REQUIREMENTS: POLYGON or MULTIPOLYGON (not point/line).
// Accepted as a FeatureCollection asset (uploaded shapefile), a single
// Feature, or a polygon drawn on the map.
var aoi = ee.FeatureCollection('projects/ee-chiararichiardi/assets/ITALYPilotSites_buffer250');

// Name of the shapefile attribute (column) to use as label for each
// AOI feature in the time-series chart. Set to null to fall back to
// system:index. Example: 'name', 'plot_id', 'sito', 'CODICE'.
var aoiLabelProperty = 'CODICE SS';

// --- Indices to compute ----------------------------------------------
// Any subset of: 'NDVI', 'NDMI'.
// e.g. ['NDVI'], ['NDMI'], or ['NDVI', 'NDMI'] for both.
var indices = ['NDVI'];

// --- Active vegetation season (day-of-year, inclusive) ----------------
// Example: 1 June - 30 September -> doyStart = 152, doyEnd = 273
var doyStart = 152;
var doyEnd   = 273;

// --- Year range -------------------------------------------------------
var startYear = 2018;
var endYear   = 2024;

// --- Quality-control thresholds --------------------------------------
var csThreshold      = 0.60; // Cloud Score+: keep pixels with cs >= this
var minImgsPerSeason = 5;    // min valid S2 scenes per pixel per season
var minValidYears    = 6;    // min valid annual composites per pixel
var pSignificance    = 0.05; // significance level for display/stats

// --- Chart y-axis range ----------------------------------------------
// Both null (default) -> Google Charts auto-scales to the data range,
// which is what you want for narrow ranges like NDVI 0.78-0.91 over
// dense vegetation.
// Set explicit numbers (e.g. chartYMin = 0, chartYMax = 1) to force the
// full theoretical index range. You can set only one of the two and
// leave the other auto-scaled.
var chartYMin = null;
var chartYMax = null;

// --- Export destination ----------------------------------------------
// Where to send the exports. One of: 'drive', 'asset', 'both'.
var exportTo = 'drive';

// Google Drive folder (root of My Drive). Created automatically if it
// does not exist. Used when exportTo includes 'drive'.
var driveFolder = 'SentinelSprings_GEE_exports';

// GEE asset folder (must already exist; create it once in the Assets
// tab). Used when exportTo includes 'asset'. Set your own username.
var assetFolder = 'users/your_username/SentinelSprings';

// Output projection (CRS) for the exported files.
// Leave null (default) to use the NATIVE Sentinel-2 UTM projection of
// the AOI, which preserves SQUARE 10 m pixels and is strongly
// recommended. Override only if a specific CRS is required by
// downstream tools, e.g. 'EPSG:3035' (LAEA Europe) for pan-European
// products. Avoid 'EPSG:4326' (rectangular pixels at mid-latitudes).
var exportCRS = null;

/* ==================== 2. AOI LOADING + DISPLAY ===================== */

var aoiFC = ee.FeatureCollection(aoi);

// Inject a stable label for each feature.
aoiFC = aoiFC.map(function (f) {
  var fallback = ee.String('feature_').cat(ee.String(f.get('system:index')));
  var lbl = aoiLabelProperty
    ? ee.Algorithms.If(f.get(aoiLabelProperty),
                       ee.String(f.get(aoiLabelProperty)), fallback)
    : fallback;
  return f.set('aoi_label', lbl);
});

var aoiGeom = aoiFC.geometry();

Map.centerObject(aoiFC);
Map.addLayer(
  aoiFC.style({color: '000000', fillColor: '00000020', width: 2}),
  {}, 'AOI features'
);
Map.addLayer(
  aoiFC.map(function (f) { return f.centroid({maxError: 1}); })
       .style({color: 'red', pointSize: 6}),
  {}, 'AOI centroids', false
);

print('--- AOI ---');
print('Number of features:', aoiFC.size());
print('Total AOI area (ha):', aoiGeom.area({maxError: 1}).divide(1e4));
print('Feature labels:', aoiFC.aggregate_array('aoi_label'));
if (aoiLabelProperty) {
  print('Labels are taken from the "' + aoiLabelProperty +
        '" attribute of the shapefile.');
} else {
  print('aoiLabelProperty is null: using system:index as label. ' +
        'Set aoiLabelProperty (e.g. "name") to use a meaningful column.');
}

/* =================== 3. SENTINEL-2 COLLECTION ====================== */

var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var CS_BAND = 'cs';

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoiGeom)
  .filterDate(ee.Date.fromYMD(startYear, 1, 1),
              ee.Date.fromYMD(endYear, 12, 31))
  .filter(ee.Filter.calendarRange(doyStart, doyEnd, 'day_of_year'))
  .linkCollection(csPlus, [CS_BAND]);

// Mask clouds/shadow with Cloud Score+, then compute the requested
// indices. Each image returned has one band per index in `indices`.
function maskAndIndices(img) {
  var masked = img.updateMask(img.select(CS_BAND).gte(csThreshold));
  var out = ee.Image().select(); // empty image, bands appended below
  if (indices.indexOf('NDVI') !== -1) {
    out = out.addBands(masked.normalizedDifference(['B8', 'B4'])
                              .rename('NDVI'));
  }
  if (indices.indexOf('NDMI') !== -1) {
    out = out.addBands(masked.normalizedDifference(['B8', 'B11'])
                              .rename('NDMI'));
  }
  return out.copyProperties(img, ['system:time_start']);
}
var s2idx = s2.map(maskAndIndices);

// Native projection of Sentinel-2 over the AOI (typically the local
// UTM zone, e.g. EPSG:32632 for Piedmont). Kept as a Projection OBJECT
// — not a derived CRS string — and imprinted on the export image via
// setDefaultProjection(). This is the robust way to preserve square
// 10 m UTM pixels in the exported GeoTIFF.
var nativeProj = s2.first().select('B8').projection();

print('--- Sentinel-2 ---');
print('Indices requested:', indices);
print('Scenes in season (all years):', s2.size());
print('Native S2 projection (used as default export CRS):', nativeProj);

/* ============== 4. ANNUAL SEASONAL p90 COMPOSITES ================== */

var years = ee.List.sequence(startYear, endYear);

var annualIdx = ee.ImageCollection.fromImages(years.map(function (y) {
  y = ee.Number(y);
  var seasonImgs = s2idx.filter(ee.Filter.calendarRange(y, y, 'year'));

  // p90 of each requested index, keeping the original band name.
  var p90 = seasonImgs.reduce(ee.Reducer.percentile([90]))
    .rename(indices); // ee.Reducer.percentile renames to <band>_p90; reset

  // Per-pixel count of valid scenes (use the first index as reference).
  var obs = seasonImgs.select(indices[0]).count().rename('obs_count');
  var validMask = obs.gte(minImgsPerSeason);

  return p90.addBands(obs)
    .updateMask(validMask)
    .set('year', y)
    .set('system:time_start', ee.Date.fromYMD(y, 7, 1).millis());
}));

// Per-pixel number of valid annual composites (same mask across indices,
// since the obs count is shared).
var validYears = annualIdx.select(indices[0]).count().rename('valid_years');

print('--- Annual composites ---');
print('Composites built:', annualIdx.size());

/* =========== 5. TIME-SERIES CHART PER AOI FEATURE ================== */
// One chart per index, one line per AOI feature.

function makeChart(indexName) {
  // y-axis viewWindow logic:
  //   - If user set chartYMin / chartYMax, use those values.
  //   - Otherwise, omit the viewWindow entirely -- Google Charts then
  //     auto-scales to the data range, which is what we want for narrow
  //     ranges (e.g. NDVI 0.78-0.91 for dense vegetation).
  var vAxis = {title: indexName + ' (p90, seasonal)'};
  if (chartYMin !== null && chartYMax !== null) {
    vAxis.viewWindow = {min: chartYMin, max: chartYMax};
  } else if (chartYMin !== null) {
    vAxis.viewWindow = {min: chartYMin};
  } else if (chartYMax !== null) {
    vAxis.viewWindow = {max: chartYMax};
  }

  var chart = ui.Chart.image.seriesByRegion({
    imageCollection: annualIdx.select(indexName),
    regions: aoiFC,
    reducer: ee.Reducer.mean(),
    scale: 10,
    seriesProperty: 'aoi_label',
    xProperty: 'system:time_start'
  })
  .setChartType('LineChart')
  .setOptions({
    title: 'Mean seasonal ' + indexName + ' (p90) per AOI feature, ' +
           startYear + '-' + endYear,
    vAxis: vAxis,
    hAxis: {title: 'Year', format: 'yyyy'},
    lineWidth: 2,
    pointSize: 5,
    interpolateNulls: false
  });
  print(chart);
}

indices.forEach(makeChart);

/* ================== 6. COVERAGE DIAGNOSTIC ========================= */

var coverageStats = validYears.reduceRegion({
  reducer: ee.Reducer.minMax()
    .combine(ee.Reducer.mean(),   '', true)
    .combine(ee.Reducer.median(), '', true),
  geometry: aoiGeom,
  scale: 10,
  maxPixels: 1e13
});
print('--- Coverage diagnostic ---');
print('valid_years per pixel (min/max/mean/median):', coverageStats);
print('Expected maximum:', (endYear - startYear + 1));

/* ============= 7. MANN-KENDALL + SEN'S SLOPE ======================= */
// Run independently per index.

function trendStats(indexName) {
  var coll = annualIdx.select(indexName);

  var afterFilter = ee.Filter.lessThan({
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });
  var joined = ee.ImageCollection(ee.Join.saveAll('after').apply({
    primary: coll, secondary: coll, condition: afterFilter
  }));

  var pairs = ee.ImageCollection(joined.map(function (current) {
    var cur = ee.Image(current);
    var curYear = ee.Number(cur.get('year'));
    var afterCol = ee.ImageCollection.fromImages(cur.get('after'));
    return afterCol.map(function (img) {
      img = ee.Image(img);
      var diff = img.subtract(cur);
      var sign = diff.gt(0).subtract(diff.lt(0)).rename('sign');
      var yearDiff = ee.Number(img.get('year')).subtract(curYear);
      var slope = diff.divide(ee.Image.constant(yearDiff)).rename('slope');
      return sign.addBands(slope);
    });
  }).flatten());

  var sStat = pairs.select('sign').reduce(ee.Reducer.sum()).rename('S');

  var n = validYears;
  var varS = n.multiply(n.subtract(1))
              .multiply(n.multiply(2).add(5))
              .divide(18)
              .rename('varS');

  var signS = sStat.gt(0).subtract(sStat.lt(0));
  var z = sStat.subtract(signS).divide(varS.sqrt()).rename('z');

  // Two-sided p-value via Zelen & Severo (1964) normal-CDF approximation
  var az = z.abs();
  var t  = ee.Image(1).divide(az.multiply(0.2316419).add(1));
  var phi = az.multiply(az).multiply(-0.5).exp().multiply(0.3989422804);
  var t2 = t.multiply(t), t3 = t2.multiply(t),
      t4 = t3.multiply(t), t5 = t4.multiply(t);
  var poly = t.multiply(0.319381530)
    .add(t2.multiply(-0.356563782))
    .add(t3.multiply(1.781477937))
    .add(t4.multiply(-1.821255978))
    .add(t5.multiply(1.330274429));
  var cdf = ee.Image(1).subtract(phi.multiply(poly));
  var pValue = ee.Image(1).subtract(cdf).multiply(2).rename('p_value');

  var nPairs = n.multiply(n.subtract(1)).divide(2);
  var tau = sStat.divide(nPairs).rename('kendall_tau');
  var sensSlope = pairs.select('slope').reduce(ee.Reducer.median())
                    .rename('sens_slope');

  // Sen's slope masked to pixels where the trend is significant
  // (p_value <= pSignificance). Non-significant pixels are masked out
  // so they appear transparent / no-data in the exported GeoTIFF.
  var sensSlopeSig = sensSlope.updateMask(pValue.lte(pSignificance))
                              .rename('sens_slope_significant');

  // Prefix band names with the index for clarity in multi-index runs.
  return sensSlope.addBands(sensSlopeSig).addBands(tau).addBands(z)
    .addBands(pValue).addBands(sStat)
    .regexpRename('^', indexName + '_');
}

// Combine results across all requested indices.
var result = indices.reduce(function (acc, idx) {
  var r = trendStats(idx);
  return acc === null ? r : acc.addBands(r);
}, null);
result = result.addBands(validYears);

/* ===================== 8. RESULTS & DISPLAY ======================== */

var enoughData = validYears.gte(minValidYears);
result = result.updateMask(enoughData).clip(aoiGeom);

// validYears also needs to be masked & clipped, otherwise it covers
// the full S2 footprint over the bounding box and looks inconsistent
// with the trend bands when both are inspected in QGIS.
var validYearsMasked = validYears.updateMask(enoughData).clip(aoiGeom);

var palette = ['a50026','d73027','f46d43','fee08b','ffffbf',
               'd9ef8b','a6d96a','66bd63','1a9850'];

indices.forEach(function (idx) {
  var slopeBand    = idx + '_sens_slope';
  var slopeSigBand = idx + '_sens_slope_significant';

  Map.addLayer(result.select(slopeBand),
    {min: -0.03, max: 0.03, palette: palette},
    idx + ' Sen slope (units/yr)', false);

  // Pre-masked band — significant pixels only.
  Map.addLayer(result.select(slopeSigBand),
    {min: -0.03, max: 0.03, palette: palette},
    idx + ' Sen slope - significant (p<' + pSignificance + ')', false);
});

Map.addLayer(validYearsMasked,
  {min: 0, max: (endYear - startYear + 1),
   palette: ['white', 'purple']},
  'Valid years per pixel', false);

// Area statistics per index.
var pxArea = ee.Image.pixelArea().divide(1e4);
indices.forEach(function (idx) {
  var slopeBand = idx + '_sens_slope';
  var pBand     = idx + '_p_value';
  var significant = result.select(pBand).lte(pSignificance);
  var slope = result.select(slopeBand);

  var areaPos = pxArea.updateMask(significant.and(slope.gt(0)));
  var areaNeg = pxArea.updateMask(significant.and(slope.lt(0)));
  var areaAll = pxArea.updateMask(slope.mask());

  print('--- ' + idx + ' trend area summary (hectares) ---');
  print(ee.Dictionary({
    'ha_with_data':            areaAll.reduceRegion({reducer: ee.Reducer.sum(),
      geometry: aoiGeom, scale: 10, maxPixels: 1e13}).get('area'),
    'ha_significant_positive': areaPos.reduceRegion({reducer: ee.Reducer.sum(),
      geometry: aoiGeom, scale: 10, maxPixels: 1e13}).get('area'),
    'ha_significant_negative': areaNeg.reduceRegion({reducer: ee.Reducer.sum(),
      geometry: aoiGeom, scale: 10, maxPixels: 1e13}).get('area')
  }));
});

/* ========================= 9. EXPORTS ============================== */
// Routes each image to Drive, Asset, or both, based on exportTo.

function exportImage(image, baseName) {
  // CRS strategy:
  //   - If exportCRS is set explicitly, pass it as a string in `crs`
  //     (the user is responsible for using a valid EPSG code).
  //   - Otherwise, IMPRINT the native S2 UTM projection directly on
  //     the image via setDefaultProjection(). This avoids GEE's
  //     "CRS could not be parsed" error that occurs when a derived
  //     ee.String CRS is passed to Export.
  var img = exportCRS ? image : image.setDefaultProjection(nativeProj);

  var base = {
    image: img,
    description: baseName,
    region: aoiGeom,
    scale: 10,
    maxPixels: 1e13
  };
  if (exportCRS) base.crs = exportCRS;

  if (exportTo === 'drive' || exportTo === 'both') {
    var pD = {};
    for (var k in base) pD[k] = base[k];
    pD.folder = driveFolder;
    pD.fileNamePrefix = baseName;
    pD.fileFormat = 'GeoTIFF';
    Export.image.toDrive(pD);
  }
  if (exportTo === 'asset' || exportTo === 'both') {
    var pA = {};
    for (var k in base) pA[k] = base[k];
    pA.assetId = assetFolder + '/' + baseName;
    Export.image.toAsset(pA);
  }
}

// 9a. Trend result — one raster PER INDEX (NDVI and NDMI are
// independent analyses on independent time series; they must NOT be
// combined into a single file).
indices.forEach(function (idx) {
  // Pick all bands for this index and drop the index prefix from the
  // band names — every file is self-contained (sens_slope,
  // sens_slope_significant, kendall_tau, z, p_value, S, valid_years).
  var indexBands = [
    idx + '_sens_slope',
    idx + '_sens_slope_significant',
    idx + '_kendall_tau',
    idx + '_z',
    idx + '_p_value',
    idx + '_S'
  ];
  var resultIdx = result.select(indexBands)
    .regexpRename('^' + idx + '_', '')
    .addBands(validYearsMasked)
    .toFloat();

  exportImage(resultIdx,
              idx + '_MK_Sen_' + startYear + '_' + endYear);
});

// 9b. Annual p90 composites, one task per year per index.
for (var y = startYear; y <= endYear; y++) {
  indices.forEach(function (idx) {
    var img = annualIdx.filter(ee.Filter.eq('year', y)).first()
                .select(idx).toFloat();
    exportImage(img, idx + '_p90_' + y);
  });
}

print('--- Exports queued ---');
print('Destination:', exportTo);
if (exportTo === 'drive' || exportTo === 'both') {
  print('Drive folder:', driveFolder);
}
if (exportTo === 'asset' || exportTo === 'both') {
  print('Asset folder:', assetFolder);
}
print('Open the Tasks panel (top-right) and click "RUN" on each task.');
