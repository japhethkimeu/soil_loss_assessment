
/** 
 * Soil Loss Estimation Using RUSLE in Google Earth Engine. https://handle.nal.usda.gov/10113/11126
 * 
 * Soil Loss Equation A = R * K * LS * C * P
 *  where, 
 *       A - Mean annual Soil Loss (metric tonnes/ha/year)
 *       R - Rainfall erosivity factor
 *       K - Soil erodibility factor
 *       LS - Slope length - steepness factor
 *       C - Land cover and management factor
 *       P - Erosion management factor
 * 
 * Author @JaphethKimeu
 * Org: FAO
*/
// ************************************************************************************************//

// import all data collections and assets
var s21 = ee.ImageCollection("COPERNICUS/S2")
var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/PENTAD")
var modis = ee.ImageCollection("MODIS/006/MCD12Q1")
var soil = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02")
var DEM = ee.Image("USGS/SRTMGL1_003")
var aoi = projects/ee-kimeu/assets/north-west_basin

// define start and end dates of analysis
var startDate = '2020-07-01';
var endDate = '2021-07-01';

// Lets compute each of the factors in the equation

// //**************** R Factor ***************// //

// filter rainfall data for the period
var filtered = chirps.filter(ee.Filter.date(startDate, endDate));
var meanAnnual = filtered.reduce(ee.Reducer.mean()).clip(aoi);

// var stats = meanAnnual.reduceRegion({
//   reducer: ee.Reducer.max(),
//   geometry: aoi,
//   scale: 5000,
//   })
// print(stats)

var palette = ['#ffffcc','#a1dab4','#41b6c4','#2c7fb8','#253494']
var visParams = {
  min:0,
  max: 10,
  palette: palette
}
// add layer to map
Map.addLayer(meanAnnual, visParams, 'Mean Precipitation')

// compute R factor (R = -8.12 + 0.562MAR, Hurni (1985). This is for Ethiopia but adapted for Somalia)
var R_factor1 = ee.Image(meanAnnual.multiply(0.562).subtract(8.12)).rename('R');

// compute R factor using
//var R_factor = ee.Image(meanAnnual.multiply(0.363).add(79)).rename('R');
// var stats = R_factor.reduceRegion({
//   reducer: ee.Reducer.min(),
//   geometry: aoi,
//   scale: 5000,
//   })
// print(stats)
// // add layer to map
Map.addLayer(R_factor1, {min: 78, max: 83, palette: palette}, 'R Factor Map', 0);

// // //**************** K Factor ***************// //

// read soil data and select the first band
soil = soil.select('b0').clip(aoi).rename('soil')    

// add layer to map
Map.addLayer(soil, {min: 0, max: 15, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'Soil', 0);

// compute K factor (Renard et al. 1997,  Wischmeier and Smith (1978))
var K_factor = soil.expression(
    "(b('soil') > 11) ? 0.0053" +
      ": (b('soil') > 10) ? 0.0170" +
        ": (b('soil') > 9) ? 0.045" +
          ": (b('soil') > 8) ? 0.050" +
            ": (b('soil') > 7) ? 0.0499" +
            ": (b('soil') > 6) ? 0.0394" +
            ": (b('soil') > 5) ? 0.0264" +
            ": (b('soil') > 4) ? 0.0423" +
            ": (b('soil') > 3) ? 0.0394" +
            ": (b('soil') > 2) ? 0.036" +
            ": (b('soil') > 1) ? 0.0341" +
            ": (b('soil') > 0) ? 0.0288" +
            ": 0")
            .rename('K').clip(aoi);    
// add layer to map
Map.addLayer(K_factor, {min: 0, max: 0.06, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'K Factor Map', 0);

// // //**************** LS Factor ***************// //

// read elevation data and compute slope
var elevation = DEM.select('elevation');
var slope_ = ee.Terrain.slope(elevation).clip(aoi);

// Convert slope to %
var slope = slope_.divide(180).multiply(Math.PI).tan().multiply(100);

// add layer to map
Map.addLayer(slope, {min: 0, max: 15, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'slope in %', 0);

// compute LS factor (LS = ((Q*M/22.13).pow(0.5) * (0.065 + 0.045S + 0.0065S*S), Wischmeier and Smith (1978))
// compute LS factor (LS = ((Q*M/100).pow(2) * (0.76 + 0.53S + 0.076S*S), Wischmeier and Smith (1978))
var ls1 = Math.sqrt(500/22.13); 
var ls2 = slope.multiply(0.53);
var ls3 = slope.pow(2).multiply(0.076);
var ls4 = ls2.add(ls3).add(0.76);
var LS_factor = ls4.multiply(ls1).rename("LS");

// add layer to map
Map.addLayer(LS_factor, {min: 3, max: 106, palette: ['a52508','ff3818','fbff18','25cdff','2f35ff','0b2dab']}, 'LS Factor Map', 0);

// //**************** C Factor ***************// //

// read sentinel 2 data and compute ndvi
s2 = s21.filterDate(startDate, endDate).sum().clip(aoi);
var ndvi_ = s2.normalizedDifference(['B8','B4'])
// add layer to map
Map.addLayer (ndvi_, {min: -1, max: 1, palette: ['FFFFFF','CC9966','CC9900', '996600', '33CC00', '009900','006600','000000']}, 'NDVI', 0);

// compute C factor (C = 0.431 - 0.805 * NDVI, De Jong (1994))
// var C_factor = ndvi_.multiply(-0.374).rename('C')

// compute C factor 
var alpha = ee.Number(-2)
var beta = ee.Number (1)

var C1 = ndvi_.multiply(alpha)
var oneImage = ee.Image(1).clip(aoi);
var C2 = oneImage.subtract(ndvi_)
var C3 = C1.divide(C2).rename('C3')
var C4 = C3.exp()

var maxC4 = C4.reduceRegion({
  geometry: aoi, 
  reducer: ee.Reducer.max(), 
  scale: 3000,
  maxPixels: 475160679
})

var C5 = maxC4.toImage().clip(aoi)
var minC4 = C4.reduceRegion({
   
  geometry: aoi, 
  reducer: ee.Reducer.min(), 
  scale: 3000,
  maxPixels: 475160679
})

var C6 = minC4.toImage().clip(aoi)
var C7 = C4.subtract(C6)
var C8 = C5.subtract(C6)

var C_factor = C7.divide(C8).rename('C')

// var stats = C_factor.reduceRegion({
//   reducer: ee.Reducer.min(),
//   geometry: aoi,
//   scale: 5000,
//   })
// print(stats)
// add layer to map
Map.addLayer (C_factor, {min: 0, max: 0.9, palette: ['FFFFFF','CC9966','CC9900', '996600', '33CC00', '009900','006600','000000']}, 'C Factor Map', 0);

// //**************** P Factor ***************// //
 
// load land cover data and filter for period of analysis
var lulc = modis.filterDate(startDate, endDate).select('LC_Type1')
        .first().clip(aoi).rename('lulc');

var lulcVis = {
  min: 1.0,
  max: 17.0,
  palette: [
    '05450a', '086a10', '54a708', '78d203', '009900', 'c6b044', 'dcd159',
    'dade48', 'fbff13', 'b6ff05', '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c',
    '69fff8', 'f9ffa4', '1c0dff'
  ],
};
// add layer to map
Map.addLayer (lulc, lulcVis, 'lulc', 0) 

// create a composite of land cover + slope
var lulc_slope = lulc.addBands(slope)

// compute P factor (Hurni et al. (2015))
var P_factor = lulc_slope.expression(
    
      "(b('lulc') < 11) ? 0.8" +
      ": (b('lulc') == 11) ? 1" +
      ": (b('lulc') == 13) ? 1" +
      ": (b('lulc') > 14) ? 1" +
      ": (b('slope') < 2) and((b('lulc')==12) or (b('lulc')==14)) ? 0.6" +
    ": (b('slope') < 5) and((b('lulc')==12) or (b('lulc')==14)) ? 0.5" +
    ": (b('slope') < 8) and((b('lulc')==12) or (b('lulc')==14)) ? 0.5" +
    ": (b('slope') < 12) and((b('lulc')==12) or (b('lulc')==14)) ? 0.6" +
    ": (b('slope') < 16) and((b('lulc')==12) or (b('lulc')==14)) ? 0.7" +
    ": (b('slope') < 20) and((b('lulc')==12) or (b('lulc')==14)) ? 0.8" +
    ": (b('slope') > 20) and((b('lulc')==12) or (b('lulc')==14)) ? 0.9" +
    ": 1"
).rename('P').clip(aoi);

// //compute P factor using Weners formula (P = 0.2 + 0.03 * S, Wener (1981))
// // var P_factor2 = slope.multiply(0.03).add(0.2).rename('P');

// // add layer to map
Map.addLayer (P_factor, {min: 0.7, max: 1, palette: ['b6ff05', '27ff87', 'c24f44', 'a5a5a5', 'ff6d4c']}, 'P Factor', 0)

// // //**************  Estimating Soil Loss ******************// //

// compute overall soil loss
var soil_loss = R_factor1.multiply(K_factor).multiply(LS_factor).multiply(C_factor).multiply(P_factor).rename("Soil Loss")
// var stats = soil_loss.reduceRegion({
//   reducer: ee.Reducer.min(),
//   geometry: aoi,
//   scale: 5000,
//   })
// print(stats)
var viz = ['490eff','12f4ff','12ff50','e5ff12','ff4812']
// add layer to map
Map.addLayer (soil_loss, {min: 1, max: 143, palette: viz}, 'Soil Loss', 0)

// create classes
var SL_class = soil_loss.expression(
    "(b('Soil Loss') < 10) ? 1" +
      ": (b('Soil Loss') < 40) ? 2" +
      ": (b('Soil Loss') < 70) ? 3" +
      ": (b('Soil Loss') < 100) ? 4" +
            ": 5")
            .rename('SL_class').clip(aoi);  
// add layer to map
Map.addLayer (SL_class, {min: 1, max: 5, palette: viz}, 'Soil Loss Class')

// compute mean soil loss
var SL_mean = soil_loss.reduceRegion({
  geometry: aoi, 
  reducer: ee.Reducer.mean(), 
  scale: 500,
  maxPixels: 475160679
})

print ("Mean Soil Loss",SL_mean.get("Soil Loss"))

// // compute area for @ class
var areaImage = ee.Image.pixelArea().addBands(
      SL_class)
var areas = areaImage.reduceRegion({
      reducer: ee.Reducer.sum().group({
      groupField: 1,
      groupName: 'class',
    }),
    geometry: aoi.geometry(),
    scale: 500,
    maxPixels: 1e10
    }); 
 
print(areas)

var classAreas = ee.List(areas.get('groups'))
 
var className = classAreas.map(function(item) {
  var areaDict = ee.Dictionary(item)
  var classNumber = ee.Number(areaDict.get('class')).format()
  return ee.List(classNumber)  
})

var Area = classAreas.map(function(item) {
  var areaDict = ee.Dictionary(item)
  var area = ee.Number(
  areaDict.get('sum')).divide(1e6).round()
  return ee.List(area) 
})

var className2 = ee.List(["Slight (<10)","Moderate (10-40)","High (40-70)","Very high (70-100)","Severe (>100)"])

print(ui.Chart.array.values(Area, 0, className2)
    .setChartType('PieChart')
    .setOptions({pointSize: 2, title: 'Soil Loss',}));


// // //********** Legend ***********// //

// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Create legend title
var legendTitle = ui.Label({
  value: 'Soil Loss (mt/ha/yr)',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette = viz;
 
// name of the legend
// var names = ["Slight","Moderate","High","Very high","Severe"];
var names = ["Slight (<10)","Moderate (10-40)","High (40-70)","Very high (70-100)","Severe (>100)"]
 
// Add color and and names
for (var i = 0; i < 5; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  
 
// add legend to map
Map.add(legend);

                //********** Links to resources used in the project *********//
// Hurni (1985) DOI 10.7892/boris.77547
// Hurni et al. (2015)
// Wischmeier and Smith (1978) https://handle.nal.usda.gov/10113/CAT79706928
// De Jong (1994)