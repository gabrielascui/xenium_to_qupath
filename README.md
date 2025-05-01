# Xenium to QuPath

This repository contains a Python script designed to facilitate the conversion of Xenium data into a geoJSON format compatible with QuPath.

## Features

- Converts Xenium output data into geoJSON formats.
- Streamlines the integration of spatial transcriptomics data into QuPath workflows.

## Requirements

- Python 3.7 or higher
- Required Python libraries (see `requirements.txt`)

## Installation

1. Clone this repository:
    ```bash
    git clone https://github.com/your-repo/xenium_to_qupath.git
    ```


## Usage

Modify the script to use the appropiated pixel size and paht the the `cells.zarr` file from the Xenium Ranger output. 

```python
# Path to the Zarr directory (adjust the path accordingly)
zarr_dir = 'cells.zarr'

# pixel size
psize = 0.2125
```


Run the script with the following command:
```bash
python xenium_to_qupath.py 
```

### Output

Output file will be a geojson file named `'exported_cells.geojson'`.

## Example

```bash
python xenium_to_qupath.py 
```


## Contact

For questions or issues, please contact [gascui@lji.org].