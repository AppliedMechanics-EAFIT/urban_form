import geopandas as gpd
from pathlib import Path

# Load the GeoJSON file
city_name = 'santafe'
PROJECT_ROOT = Path(__file__).resolve().parents[1]
input_path = PROJECT_ROOT / "labeled_data" / "original" / f"urban_patterns_clustered_{city_name}.geojson"
output_path = PROJECT_ROOT / "labeled_data" / "with_centroids" / f"urban_patterns_clustered_{city_name}.geojson"


gdf = gpd.read_file(input_path)

# Ensure the original CRS is WGS84
if gdf.crs is None:
    gdf.set_crs(epsg=4326, inplace=True)

# Step 1: Reproject to a suitable projected CRS
# Option A: Automatically estimate the best UTM zone
gdf_projected = gdf.to_crs(gdf.estimate_utm_crs())

# Step 2: Compute centroids in the projected CRS
centroids_projected = gdf_projected.geometry.centroid

# Step 3: Transform centroids back to WGS84
centroids_wgs84 = gpd.GeoSeries(centroids_projected, crs=gdf_projected.crs).to_crs(epsg=4326)

# Step 4: Add centroid coordinates to properties
gdf["center_lon"] = centroids_wgs84.x
gdf["center_lat"] = centroids_wgs84.y

# Optional: Store as a GeoJSON Point
gdf["center"] = gdf.apply(
    lambda row: {
        "type": "Point",
        "coordinates": [row["center_lat"], row["center_lon"]],
    },
    axis=1
)

# Save the updated GeoJSON
output_path.parent.mkdir(parents=True, exist_ok=True)
gdf.to_file(output_path, driver="GeoJSON")

print("Centroids successfully added!")
