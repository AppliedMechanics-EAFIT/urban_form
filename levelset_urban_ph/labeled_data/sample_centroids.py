import json
import random
from pathlib import Path

PATTERNS = ["cul_de_sac", "gridiron", "organico", "hibrido"]
N = 30000
random.seed(10)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
CENTROIDS_DIR = PROJECT_ROOT / "labeled_data" / "with_centroids"

results = []
id_counters = {}

for filepath in sorted(CENTROIDS_DIR.glob("urban_patterns_clustered_*.geojson")):
    city = filepath.stem.split("clustered_", 1)[1]
    with filepath.open() as f:
        data = json.load(f)

    by_pattern = {p: [] for p in PATTERNS}
    for feat in data["features"]:
        p = feat["properties"].get("original_pattern")
        if p in by_pattern:
            by_pattern[p].append(feat)

    for pattern, feats in by_pattern.items():
        sample = random.sample(feats, min(N, len(feats)))

        for feat in sample:
            props = feat["properties"]

            key = (city, pattern)
            id_counters[key] = id_counters.get(key, 0) + 1
            sample_id = f"{city}_{pattern}_{id_counters[key]}"

            results.append({
                "id": sample_id,
                "city": city,
                "coordinates": [props["center_lat"], props["center_lon"]],
                "original_pattern": pattern,
                "geometry": feat["geometry"],
            })

output_path = CENTROIDS_DIR / "sample_everything.json"
with output_path.open("w") as f:
    json.dump(results, f, indent=2)

print(f"Saved {len(results)} records to {output_path}")
