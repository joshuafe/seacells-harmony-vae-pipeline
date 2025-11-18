
#!/bin/bash

set -e

# Create output directory
mkdir -p single_sample_analysis_R_output

# Process normal samples
for fcs_file in normalsgated/*.fcs; do
    sample_name=$(basename "$fcs_file" .fcs)
    echo "Processing normal sample: $sample_name"
    output_dir="single_sample_analysis_R_output/normal_$sample_name"
    Rscript single_sample_analysis.R "$fcs_file" "$output_dir"
done

# Process abnormal samples
for fcs_file in abnormalfcs_gated/*.fcs; do
    sample_name=$(basename "$fcs_file" .fcs)
    echo "Processing abnormal sample: $sample_name"
    output_dir="single_sample_analysis_R_output/abnormal_$sample_name"
    Rscript single_sample_analysis.R "$fcs_file" "$output_dir"
done

echo "All samples processed."
