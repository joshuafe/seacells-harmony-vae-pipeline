"""
Inspects an FCS file and prints its marker channel names (v2).
This version is more robust and inspects the data columns directly.
"""
import sys
from fcsparser import parse
from pathlib import Path
import re

def inspect_fcs_markers_v2(fcs_path):
    try:
        fcs_path = Path(fcs_path)
        if not fcs_path.exists():
            print(f"Error: File not found at {fcs_path}")
            return

        print(f"Inspecting file: {fcs_path.name}")
        
        # Use reformat_meta=False to see the raw metadata keys
        meta, data = parse(str(fcs_path), reformat_meta=False)
        
        print("\n--- Raw Metadata Keys ---")
        print(list(meta.keys()))

        print("\n--- Data Columns ---")
        print(data.columns)

        # The column names are usually the most reliable source of markers
        # They might be in the format 'Comp-BV421-A :: CD16+TIGIT'
        # Or just 'CD16+TIGIT'
        
        final_markers = []
        print("\n" + "="*40)
        print("PARSING MARKERS FROM DATA COLUMNS")
        print("="*40)

        for col in data.columns:
            # Try to extract marker from 'Fluorochrome :: Marker' format
            parts = col.split('::')
            if len(parts) > 1:
                marker_name = parts[-1].strip()
            else:
                marker_name = col.strip()

            # Filter out common non-marker channels
            if re.search(r'FSC|SSC|Time|Viability', marker_name, re.IGNORECASE):
                print(f"Ignoring non-marker channel: '{col}' -> '{marker_name}'")
                continue
            
            print(f"Found marker: '{col}' -> '{marker_name}'")
            final_markers.append(marker_name)
            
        print("\n" + "="*40)
        print("FINAL PARSED MARKER LIST")
        print("="*40)
        print(f"Total markers found: {len(final_markers)}")
        # Use repr() to make it easy to copy-paste into a Python list
        print(repr(final_markers))

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_fcs_v2.py <path_to_fcs_file>")
    else:
        inspect_fcs_markers_v2(sys.argv[1])
