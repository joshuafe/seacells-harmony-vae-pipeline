import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

def analyze_disambiguation_profiles(csv_path, output_dir):
    """
    Loads disambiguated archetype profiles, generates a heatmap,
    and prints a biological interpretation of each archetype.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    try:
        df = pd.read_csv(csv_path, index_col=0)
    except FileNotFoundError:
        print(f"Error: File not found at {csv_path}")
        return

    print("="*80)
    print("BIOLOGICAL INTERPRETATION OF DISAMBIGUATED ARCHETYPES")
    print("="*80)

    # Generate and save heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(df, cmap='RdBu_r', center=0, xticklabels=True, yticklabels=True)
    plt.title('Disambiguated Archetype Marker Profiles')
    plt.xlabel('Un-merged Markers')
    plt.ylabel('Archetypes')
    heatmap_path = output_dir / "disambiguated_profiles_heatmap.png"
    plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved heatmap to {heatmap_path}\n")

    # Interpret each archetype
    for index, row in df.iterrows():
        # Get top 5 positive and top 2 negative markers
        top_markers = row.sort_values(ascending=False)
        top_pos = top_markers.head(5)
        top_neg = top_markers.tail(2)

        print(f"\n--- {index} ---")
        
        # Basic interpretation logic
        is_t_cell = top_pos.get('CD3', -1) > 0.5
        is_b_cell = top_pos.get('CD19', -1) > 0.5
        is_nk_cell = top_pos.get('CD56', -1) > 0.5 or top_pos.get('CD16', -1) > 0.5
        
        is_cd4 = top_pos.get('CD4', -1) > 0.5
        is_cd8 = top_pos.get('CD8', -1) > 0.5
        
        is_naive = top_pos.get('CD45RA', -1) > 0.5 and top_pos.get('CD197', -1) > 0.5
        is_memory = top_pos.get('CD45RO', -1) > 0.5
        is_em = is_memory and top_markers.get('CD197', 1) < 0.0
        is_cm = is_memory and top_markers.get('CD197', -1) > 0.0

        interpretation = "Uncertain"
        if is_t_cell:
            if is_cd4:
                if is_naive: interpretation = "Naive CD4+ T-cell"
                elif is_cm: interpretation = "Central Memory CD4+ T-cell"
                elif is_em: interpretation = "Effector Memory CD4+ T-cell"
                else: interpretation = "CD4+ T-cell"
            elif is_cd8:
                if is_naive: interpretation = "Naive CD8+ T-cell"
                elif is_cm: interpretation = "Central Memory CD8+ T-cell"
                elif is_em: interpretation = "Effector Memory CD8+ T-cell"
                else: interpretation = "CD8+ T-cell"
            else:
                interpretation = "T-cell (CD4/CD8 low)"
        elif is_b_cell:
            interpretation = "B-cell"
        elif is_nk_cell:
            interpretation = "NK cell"
        
        print(f"Interpretation: {interpretation}")
        print("Top Positive Markers:")
        for marker, value in top_pos.items():
            print(f"  - {marker}: {value:.2f}")
        print("Top Negative Markers:")
        for marker, value in top_neg.items():
            print(f"  - {marker}: {value:.2f}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        csv_path = sys.argv[1]
        output_dir = Path(csv_path).parent
        analyze_disambiguation_profiles(csv_path, output_dir)
    else:
        print("Usage: python analyze_disambiguation.py <path_to_csv>")
