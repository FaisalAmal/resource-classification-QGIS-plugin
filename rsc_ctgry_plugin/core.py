import geopandas as gpd
from shapely.geometry import MultiPoint, Polygon, MultiPolygon, Point
from shapely.ops import unary_union
from scipy.spatial import Delaunay
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import os

def calculate_nearest_neighbor_distances(points_gdf):
    """Calculate nearest neighbor distance for each drillhole"""
    from scipy.spatial import cKDTree
    
    coords = np.array(list(zip(points_gdf.geometry.x, points_gdf.geometry.y)))
    
    if len(coords) < 2:
        return np.zeros(len(coords))

    tree = cKDTree(coords)
    distances, indices = tree.query(coords, k=2)
    return distances[:, 1]

def calculate_triangle_area(p1, p2, p3):
    """Calculate area of a triangle given three points"""
    v1 = p2 - p1
    v2 = p3 - p1
    area = 0.5 * abs(np.cross(v1, v2))
    return area

def is_valid_triangle(p1, p2, p3, max_area, max_edge, edge_tolerance=1.5):
    """Check if triangle is valid based on area and edge constraints"""
    area = calculate_triangle_area(p1, p2, p3)
    if area > max_area:
        return False
    
    edge1 = np.linalg.norm(p2 - p1)
    edge2 = np.linalg.norm(p3 - p2)
    edge3 = np.linalg.norm(p1 - p3)
    edges = sorted([edge1, edge2, edge3])
    
    # At least 2 edges within limit, longest can be diagonal
    if edges[0] <= max_edge and edges[1] <= max_edge:
        if edges[2] <= max_edge * edge_tolerance:
            return True
    
    return False

def find_triangle_based_clusters(points_gdf, max_distance, max_area, min_drillholes=4, edge_tolerance=1.5):
    """Find clusters using triangle connectivity"""
    coords = np.array(list(zip(points_gdf.geometry.x, points_gdf.geometry.y)))
    n_points = len(coords)
    
    if n_points < min_drillholes:
        return np.full(n_points, -1)
    
    tri = Delaunay(coords)
    
    # Validate triangles
    valid_triangle_indices = []
    
    for simplex in tri.simplices:
        p1, p2, p3 = coords[simplex]
        
        area = calculate_triangle_area(p1, p2, p3)
        if area > max_area:
            continue
        
        edge1 = np.linalg.norm(p2 - p1)
        edge2 = np.linalg.norm(p3 - p2)
        edge3 = np.linalg.norm(p1 - p3)
        edges = sorted([edge1, edge2, edge3])
        
        if edges[0] <= max_distance and edges[1] <= max_distance:
            if edges[2] <= max_distance * edge_tolerance:
                valid_triangle_indices.append(simplex)
    
    if not valid_triangle_indices:
        return np.full(n_points, -1)
    
    # Build adjacency from valid triangles
    adjacency = np.zeros((n_points, n_points), dtype=int)
    for simplex in valid_triangle_indices:
        adjacency[simplex[0], simplex[1]] = 1
        adjacency[simplex[1], simplex[0]] = 1
        adjacency[simplex[0], simplex[2]] = 1
        adjacency[simplex[2], simplex[0]] = 1
        adjacency[simplex[1], simplex[2]] = 1
        adjacency[simplex[2], simplex[1]] = 1
    
    # Find connected components
    graph = csr_matrix(adjacency)
    n_components, labels = connected_components(csgraph=graph, directed=False)
    
    # Filter by minimum size
    valid_labels = labels.copy()
    for cluster_id in range(n_components):
        if np.sum(labels == cluster_id) < min_drillholes:
            valid_labels[labels == cluster_id] = -1
    
    # Renumber
    unique_valid = np.unique(valid_labels[valid_labels >= 0])
    label_map = {old: new for new, old in enumerate(unique_valid)}
    label_map[-1] = -1
    final_labels = np.array([label_map[label] for label in valid_labels])
    
    return final_labels

def calculate_adaptive_alpha(max_distance, area_factor, edge_tolerance):
    """
    Calculate alpha radius that adapts to triangle constraints.
    """
    # Estimate circumradius for triangle at constraint limits
    max_edge = max_distance * edge_tolerance
    
    # Circumradius is roughly proportional to edge length
    # Scale by area_factor to respect area constraint
    alpha = max_edge * 0.6 * np.sqrt(area_factor)
    
    return alpha

def alpha_shape(points, alpha_radius):
    """
    Create alpha shape using circumradius constraint.
    """
    if len(points) < 3:
        return MultiPoint(points).convex_hull
    
    tri = Delaunay(points)
    
    def circumradius(p1, p2, p3):
        """Calculate circumradius of triangle"""
        a = np.linalg.norm(p2 - p1)
        b = np.linalg.norm(p3 - p2)
        c = np.linalg.norm(p1 - p3)
        area = 0.5 * abs(np.cross(p2 - p1, p3 - p1))
        if area < 1e-10:
            return float('inf')
        return (a * b * c) / (4.0 * area)
    
    # Filter triangles by circumradius
    valid_triangles = []
    for simplex in tri.simplices:
        p1, p2, p3 = points[simplex]
        r = circumradius(p1, p2, p3)
        
        # Keep triangles with circumradius <= alpha
        if r <= alpha_radius:
            triangle = Polygon([p1, p2, p3])
            valid_triangles.append(triangle)
    
    if not valid_triangles:
        # If alpha too restrictive, use convex hull
        return MultiPoint(points).convex_hull
    
    # Union all valid triangles
    alpha_poly = unary_union(valid_triangles)
    return alpha_poly

def classify_drillholes(drillholes, measured_distance=71, indicated_distance=121, 
                        min_drillholes=4, area_factor=0.6, edge_tolerance=1.5):
    """Classify drillholes using triangle-based clustering"""
    
    nn_distances = calculate_nearest_neighbor_distances(drillholes)
    drillholes['nn_distance'] = nn_distances
    
    # Initial classification
    drillholes['category'] = 'Inferred'
    drillholes.loc[nn_distances < indicated_distance, 'category'] = 'Indicated'
    drillholes.loc[nn_distances < measured_distance, 'category'] = 'Measured'
    
    # Initialize cluster columns for ALL drillholes
    drillholes['measured_cluster'] = -1
    drillholes['indicated_cluster'] = -1
    
    # Triangle clustering for Measured
    measured_mask = drillholes['category'] == 'Measured'
    
    if measured_mask.sum() >= min_drillholes:
        max_area = measured_distance ** 2 * area_factor
        
        measured_subset = drillholes[measured_mask].copy()
        measured_clusters = find_triangle_based_clusters(
            measured_subset, measured_distance, max_area, 
            min_drillholes, edge_tolerance
        )
        
        # Assign clusters back to original dataframe using index
        drillholes.loc[measured_mask, 'measured_cluster'] = measured_clusters
    
    # Triangle clustering for Indicated (includes Measured holes)
    indicated_mask = drillholes['category'].isin(['Measured', 'Indicated'])
    
    if indicated_mask.sum() >= min_drillholes:
        max_area = indicated_distance ** 2 * area_factor
        
        indicated_subset = drillholes[indicated_mask].copy()
        indicated_clusters = find_triangle_based_clusters(
            indicated_subset, indicated_distance, max_area, 
            min_drillholes, edge_tolerance
        )
        
        # Assign clusters back to original dataframe using index
        drillholes.loc[indicated_mask, 'indicated_cluster'] = indicated_clusters
    
    # Downgrade invalid clusters
    drillholes.loc[(drillholes['category'] == 'Measured') & (drillholes['measured_cluster'] == -1), 'category'] = 'Indicated'
    drillholes.loc[(drillholes['category'] == 'Indicated') & (drillholes['indicated_cluster'] == -1), 'category'] = 'Inferred'
    
    # Re-cluster Indicated after downgrade
    indicated_mask = drillholes['category'].isin(['Measured', 'Indicated'])
    if indicated_mask.sum() >= min_drillholes:
        max_area = indicated_distance ** 2 * area_factor
        indicated_subset = drillholes[indicated_mask].copy()
        indicated_clusters = find_triangle_based_clusters(
            indicated_subset, indicated_distance, max_area, 
            min_drillholes, edge_tolerance
        )
        drillholes.loc[indicated_mask, 'indicated_cluster'] = indicated_clusters
    
    drillholes.loc[(drillholes['category'] == 'Indicated') & (drillholes['indicated_cluster'] == -1), 'category'] = 'Inferred'
    
    return drillholes

def generate_alpha_shape_polygons(drillholes, measured_distance, indicated_distance,
                                  area_factor, edge_tolerance,
                                  measured_alpha=None, indicated_alpha=None,
                                  measured_buffer=25, indicated_buffer=50):
    """
    Generate polygons using alpha shapes that adapt to triangle constraints.
    If alpha not specified, automatically calculated from triangle parameters.
    """
    
    all_polygons = []
    
    # Process Measured
    measured_valid = drillholes[drillholes['measured_cluster'] >= 0].copy()
    
    if len(measured_valid) > 0:
        # Auto-calculate alpha from triangle constraints if not specified
        if measured_alpha is None:
            measured_alpha_val = calculate_adaptive_alpha(measured_distance, area_factor, edge_tolerance)
        else:
            measured_alpha_val = measured_alpha
        
        for cluster_id in measured_valid['measured_cluster'].unique():
            if cluster_id < 0:
                continue
                
            cluster_holes = measured_valid[measured_valid['measured_cluster'] == cluster_id]
            coords = np.array(list(zip(cluster_holes.geometry.x, cluster_holes.geometry.y)))
            
            # Alpha shape generates the polygon
            alpha_poly = alpha_shape(coords, measured_alpha_val)
            buffered_poly = alpha_poly.buffer(measured_buffer)
            
            all_polygons.append({
                'geometry': buffered_poly,
                'category': 'Measured',
                'cluster_id': int(cluster_id),
                'drillholes': len(cluster_holes)
            })
    
    # Process Indicated
    indicated_valid = drillholes[drillholes['indicated_cluster'] >= 0].copy()
    
    if len(indicated_valid) > 0:
        # Auto-calculate alpha from triangle constraints if not specified
        if indicated_alpha is None:
            indicated_alpha_val = calculate_adaptive_alpha(indicated_distance, area_factor, edge_tolerance)
        else:
            indicated_alpha_val = indicated_alpha
        
        measured_geoms = [p['geometry'] for p in all_polygons if p['category'] == 'Measured']
        measured_union = unary_union(measured_geoms) if measured_geoms else None
        
        for cluster_id in indicated_valid['indicated_cluster'].unique():
            if cluster_id < 0:
                continue
                
            cluster_holes = indicated_valid[indicated_valid['indicated_cluster'] == cluster_id]
            coords = np.array(list(zip(cluster_holes.geometry.x, cluster_holes.geometry.y)))
            
            # Alpha shape generates the polygon
            alpha_poly = alpha_shape(coords, indicated_alpha_val)
            buffered_poly = alpha_poly.buffer(indicated_buffer)
            
            if measured_union:
                buffered_poly = buffered_poly.difference(measured_union)
            
            if not buffered_poly.is_empty:
                all_polygons.append({
                    'geometry': buffered_poly,
                    'category': 'Indicated',
                    'cluster_id': int(cluster_id),
                    'drillholes': len(cluster_holes)
                })
    
    if not all_polygons:
        return gpd.GeoDataFrame()
    
    result_gdf = gpd.GeoDataFrame(all_polygons, crs=drillholes.crs)
    result_gdf['area_m2'] = result_gdf.geometry.area
    result_gdf['area_ha'] = result_gdf['area_m2'] / 10000
    
    return result_gdf

def generate_resource_polygons(input_data, output_path=None,
                               measured_distance=71, indicated_distance=121,
                               measured_alpha=None, indicated_alpha=None,
                               measured_buffer=25, indicated_buffer=50,
                               min_drillholes=4, area_factor=0.6, edge_tolerance=1.5):
    """
    Generate resource polygons.
    input_data can be a filepath (str) or a GeoDataFrame.
    """
    try:
        if isinstance(input_data, str):
            drillholes = gpd.read_file(input_data)
        elif isinstance(input_data, gpd.GeoDataFrame):
            drillholes = input_data.copy()
        else:
            raise ValueError("Input data must be a file path string or a GeoDataFrame")
            
        if drillholes.empty:
            return gpd.GeoDataFrame(), drillholes

        # Step 1: Triangle-based clustering
        classified = classify_drillholes(
            drillholes, measured_distance, indicated_distance,
            min_drillholes, area_factor, edge_tolerance
        )
        
        # Step 2: Alpha shape polygon generation
        result_gdf = generate_alpha_shape_polygons(
            classified, measured_distance, indicated_distance,
            area_factor, edge_tolerance,
            measured_alpha, indicated_alpha,
            measured_buffer, indicated_buffer
        )
        
        # Save if output path provided
        if output_path and not result_gdf.empty:
            # result_gdf.to_file(output_path, driver='GPKG')
            # For QGIS plugins, sometimes it's better to just return the GDF
            # and let the plugin handle loading/saving, but the original code saved it.
            # We will keep saving if output_path is present.
            result_gdf.to_file(output_path, driver='GPKG')
            
            # Optional: Save classified points
            # drill_out = output_path.replace('.gpkg', '_drillholes.gpkg')
            # classified.to_file(drill_out, driver='GPKG')

        return result_gdf, classified
        
    except Exception as e:
        # In a real plugin it's good to log this
        print(f"Core processing error: {e}")
        raise e
