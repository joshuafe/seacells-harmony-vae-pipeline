import pandas as pd
import sys

if len(sys.argv) < 2:
    print("Usage: python show_csv.py <path_to_csv>")
    sys.exit(1)

csv_path = sys.argv[1]

try:
    df = pd.read_csv(csv_path)
    print(df.to_string())
except FileNotFoundError:
    print(f"Error: File not found at {csv_path}")
except Exception as e:
    print(f"An error occurred: {e}")
