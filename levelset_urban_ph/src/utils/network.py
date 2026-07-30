"""Street network extraction and rasterization."""
import numpy as np
import osmnx as ox
from skimage.draw import line as draw_line


def fetch_network(center_point, dist_m):
    G = ox.graph_from_point(center_point, dist=dist_m, network_type='drive', retain_all=True)
    G_proj = ox.project_graph(G)
    nodes, edges = ox.graph_to_gdfs(G_proj)
    x_min = nodes.geometry.x.min()
    y_min = nodes.geometry.y.min()
    x_max = nodes.geometry.x.max()
    y_max = nodes.geometry.y.max()
    print(f'  Nodes: {len(nodes)}, Edges: {len(edges)}')
    print(f'  Bounding box (m): Δx={x_max-x_min:.0f}, Δy={y_max-y_min:.0f}')
    return nodes, edges, x_min, x_max, y_min, y_max


def rasterize_streets(edges, nodes, img_size, x_min, x_max, y_min, y_max, line_width=2):
    canvas = np.ones((img_size, img_size), dtype=np.float32)

    def to_pixel(x, y):
        px = int((x - x_min) / (x_max - x_min) * (img_size - 1))
        py = int((y_max - y) / (y_max - y_min) * (img_size - 1))
        return np.clip(px, 0, img_size - 1), np.clip(py, 0, img_size - 1)

    for _, edge in edges.iterrows():
        geom = edge.geometry
        if geom is None:
            continue
        coords = list(geom.coords)
        for i in range(len(coords) - 1):
            x0, y0 = to_pixel(*coords[i][:2])
            x1, y1 = to_pixel(*coords[i + 1][:2])
            rr, cc = draw_line(y0, x0, y1, x1)
            for dy in range(-line_width, line_width + 1):
                for dx in range(-line_width, line_width + 1):
                    nr = np.clip(rr + dy, 0, img_size - 1)
                    nc = np.clip(cc + dx, 0, img_size - 1)
                    canvas[nr, nc] = 0.0

    print(f'Rasterized: {canvas.shape}, street pixels: {int((canvas == 0).sum())}')
    return canvas
