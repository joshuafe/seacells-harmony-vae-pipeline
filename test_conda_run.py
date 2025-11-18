
import sys
import os

print("Hello from test_conda_run.py to stdout")
print("Hello from test_conda_run.py to stderr", file=sys.stderr)

with open("test_conda_run_output.txt", "w") as f:
    f.write("This is a test output file.\n")

print("Script finished successfully.")
