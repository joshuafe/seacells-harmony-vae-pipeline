"""
Inspects an FCS file and prints its marker channel names.
"""
import sys
from fcsparser import parse
from pathlib import Path

def inspect_fcs_markers(fcs_path):
    try:
        fcs_path = Path(fcs_path)
        if not fcs_path.exists():
            print(f"Error: File not found at {fcs_path}")
            return

        print(f"Inspecting file: {fcs_path.name}")
        
        meta, data = parse(str(fcs_path), reformat_meta=True)
        
        # Channel names are often in the format 'Comp-Fluorochrome-A :: Marker'
        # We want to extract the 'Marker' part.
        # The 'channels' metadata is the most reliable source.
        channel_info = meta['channels']
        
        print("\n" + "="*40)
        print("MARKER PANEL FROM FCS FILE METADATA")
        print("="*40)
        
        markers = []
        for i in range(len(channel_info['$P_NAME'])):
            p_name = channel_info[f'$P{i+1}N']
            p_stain = channel_info.get(f'$P{i+1}S', 'N/A')
            
            # The stain description is usually the most reliable marker name
            if p_stain != 'N/A' and p_stain.strip() != '':
                marker_name = p_stain
            else:
                marker_name = p_name
            
            # Clean up common prefixes/suffixes
            marker_name = marker_name.replace('Comp-', '').replace('-A', '')
            
            print(f"Channel {i+1}: {p_name:<25} -> Stain: '{p_stain}' -> Parsed: '{marker_name}'")
            markers.append(marker_name)
            
        print("\n" + "="*40)
        print("FINAL PARSED MARKER LIST")
        print("="*40)
        print(f"Total markers found: {len(markers)}")
        # Use repr() to make it easy to copy-paste into a Python list
        print(repr(markers))

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python inspect_fcs.py <path_to_fcs_file>")
    else:
        inspect_fcs_markers(sys.argv[1])
