#!/usr/bin/env python3

import zarr
import json
import numpy as np
import uuid  # For generating unique IDs

# Path to the Zarr directory (adjust the path accordingly)
zarr_dir = 'cells.zarr'

# pixel size
psize = 0.2125

# Open the Zarr root group
z = zarr.open(zarr_dir, mode='r')

# Load the `cell_id` dataset
cell_id_data = z['cell_id'][:]  # Shape (194412, 2)

# Function to convert NumPy types to native Python types
def convert_to_native(obj):
    """Convert any NumPy types to native Python types."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert arrays to lists
    elif isinstance(obj, (np.int32, np.int64, np.uint32, np.uint64)):
        return int(obj)  # Convert NumPy integers to Python int
    elif isinstance(obj, (np.float32, np.float64)):
        return float(obj)  # Convert NumPy floats to Python float
    elif isinstance(obj, dict):
        return {k: convert_to_native(v) for k, v in obj.items()}  # Recursively apply for dict
    elif isinstance(obj, list):
        return [convert_to_native(i) for i in obj]  # Recursively apply for lists
    return obj


#helper function to convert hex to Xenium string
# https://www.10xgenomics.com/support/software/xenium-onboard-analysis/1.9/analysis/xoa-output-zarr#cellID
def shiftCharacters(c):
    return chr(ord('a')+int(c,16))

# Helper function to process each polygon set and link it to cell_id
def process_polygon_set(set_index):
    features_by_cell_id = {}
    
    # Access specific arrays inside the current polygon set group
    vertices = z[f'polygon_sets/{set_index}/vertices'][:]  # Vertices for polygon set
    num_vertices = z[f'polygon_sets/{set_index}/num_vertices'][:]  # Number of vertices for each polygon
    cell_index = z[f'polygon_sets/{set_index}/cell_index'][:]  # Cell index for each polygon

    # Iterate through the polygons and convert them to GeoJSON features
    for i in range(len(vertices)):
        n_vertices = int(num_vertices[i])  # Ensure this is a native int
        
        if n_vertices > 0:  # Ensure there's a valid polygon
            # Extract the valid vertices (pairs of x, y)
            vertex_data = vertices[i][:n_vertices * 2]  # Take only the valid part of the array
            valid_vertices = [[float(vertex_data[j]/psize), float(vertex_data[j + 1]/psize)] for j in range(0, len(vertex_data), 2)]
            
            # GeoJSON polygons need to be closed (last point = first point)
            if valid_vertices[0] != valid_vertices[-1]:
                valid_vertices.append(valid_vertices[0])  # Close the polygon
            
            # Get the corresponding cell_id using the cell_index
            cell_id = int(cell_id_data[cell_index[i]][0])  # First value in each pair is the unique cell_id
            cellhex = format(cell_id,'08x')
            cellIDshifted = ''.join(shiftCharacters(c) for c in cellhex)
            cellIDfinal = cellIDshifted+'-'+str(cell_id_data[cell_index[i]][1])
            
            # Check if the cell_id already exists in the feature collection
            if cell_id not in features_by_cell_id:
                # Create a new feature for this cell_id
                features_by_cell_id[cell_id] = {
                    "type": "Feature",
                    "id": str(uuid.uuid4()),  # Generate a unique ID for each feature
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": []  # Initialize coordinates for cell
                    },
                    "nucleusGeometry": {
                        "type": "Polygon",
                        "coordinates": []  # Initialize coordinates for nucleus
                    },
                    "properties": {
                        "cell_id": cell_id,  # Only the cell_id is kept in the properties
                        "objectType": "cell",
                        "name": cellIDfinal
                    }
                }
            
            # Add the polygon to the correct geometry (nucleus for set 0, cell for set 1)
            if set_index == 0:
                features_by_cell_id[cell_id]["nucleusGeometry"]["coordinates"].append(valid_vertices)
            else:
                features_by_cell_id[cell_id]["geometry"]["coordinates"].append(valid_vertices)

    return features_by_cell_id

# Process both polygon sets (0 for nucleus, 1 for cell)
features_set_0 = process_polygon_set(0)
features_set_1 = process_polygon_set(1)

# Merge features from both sets based on cell_id
all_features = {}

# Merge nucleus polygons (set 0)
for cell_id, feature in features_set_0.items():
    all_features[cell_id] = feature

# Merge cell polygons (set 1)
for cell_id, feature in features_set_1.items():
    if cell_id in all_features:
        # If the cell_id exists, add the cell boundary geometry
        all_features[cell_id]["geometry"]["coordinates"].extend(feature["geometry"]["coordinates"])
    else:
        # If the cell_id doesn't exist (which is unlikely), create a new feature
        all_features[cell_id] = feature

# Convert the merged features into a list for FeatureCollection
features_list = list(all_features.values())

# Create a FeatureCollection to store all features
feature_collection = {
    "type": "FeatureCollection",
    "features": features_list
}

# Ensure all data is converted to native types for JSON serialization
native_feature_collection = convert_to_native(feature_collection)

# Write the FeatureCollection to a GeoJSON file
with open('exported_cells.geojson', 'w') as geojson_file:
    json.dump(native_feature_collection, geojson_file, separators=(',', ':'), indent=4)

print("GeoJSON file with separate nucleus and cell boundaries created successfully.")
