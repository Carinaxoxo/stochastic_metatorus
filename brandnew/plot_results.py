import pandas as pd
import matplotlib.pyplot as plt
import re

def parse_results_file(filepath):
    """
    Parses the formatted experiment_results.txt file to extract data tables.
    """
    try:
        with open(filepath, 'r') as f:
            text = f.read()
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        print("Please run your C++ experiment first to generate the results file.")
        return None

    # Split the file content by the main experiment header
    experiments = text.strip().split('--- Starting Experiment for Torus')

    all_data = {}

    for experiment in experiments:
        if not experiment.strip():
            continue

        # Find the n and k values for this experiment block
        nk_match = re.search(r'n=(\d+), k=(\d+)', experiment)
        if not nk_match:
            continue
        n, k = nk_match.groups()
        config_key = f"n={n}, k={k}"

        # Find the final results table within the block
        results_table_match = re.search(r'--- Final Results for n=\d+, k=\d+ ---(.*?)---', experiment, re.DOTALL)
        if not results_table_match:
            continue

        table_text = results_table_match.group(1).strip()

        lines = table_text.split('\n')
        data_rows = []

        # Start parsing after the header lines
        for line in lines[3:]:
            # Use regex to find all floating point or integer numbers in the line
            numbers = re.findall(r'[\d\.]+', line)
            if len(numbers) == 9: # Expecting 9 numbers per data row
                # Convert string numbers to float
                data_rows.append([float(num) for num in numbers])

        # Create a pandas DataFrame
        df = pd.DataFrame(data_rows, columns=[
            'FaultRatio', 'AvgBfsLen',
            'Brute_Success', 'Brute_AvgLen',
            'Directed_Success', 'Directed_AvgLen',
            'Strat2_Success', 'Strat2_AvgLen',
            'Strat3_Success', 'Strat3_AvgLen'
        ])
        all_data[config_key] = df

    return all_data

def create_plots(config_key, data):
    """
    Generates and saves the two graphs for a given dataset.
    """
    # --- Plot 1: Success Rate ---
    plt.figure(figsize=(12, 7))
    plt.plot(data['FaultRatio'], data['Brute_Success'], marker='o', label='Brute')
    plt.plot(data['FaultRatio'], data['Directed_Success'], marker='s', label='Directed (Simple)')
    plt.plot(data['FaultRatio'], data['Strat2_Success'], marker='^', label='Strategic (2-Class)')
    plt.plot(data['FaultRatio'], data['Strat3_Success'], marker='x', label='Strategic (3-Class)')

    plt.title(f'Algorithm Success Rate vs. Fault Ratio ({config_key})')
    plt.xlabel('Fault Ratio (%)')
    plt.ylabel('Success Rate (%)')
    plt.grid(True)
    plt.legend()
    plt.ylim(0, 101) # Set y-axis from 0 to 100%
    plt.savefig(f"success_rate_{config_key.replace(', ', '_')}.png")
    plt.show()

    # --- Plot 2: Average Path Length ---
    plt.figure(figsize=(12, 7))
    plt.plot(data['FaultRatio'], data['AvgBfsLen'], marker='d', linestyle='--', color='black', label='Optimal (BFS)')
    plt.plot(data['FaultRatio'], data['Brute_AvgLen'], marker='o', label='Brute')
    plt.plot(data['FaultRatio'], data['Directed_AvgLen'], marker='s', label='Directed (Simple)')
    plt.plot(data['FaultRatio'], data['Strat2_AvgLen'], marker='^', label='Strategic (2-Class)')
    plt.plot(data['FaultRatio'], data['Strat3_AvgLen'], marker='x', label='Strategic (3-Class)')

    plt.title(f'Average Path Length vs. Fault Ratio ({config_key})')
    plt.xlabel('Fault Ratio (%)')
    plt.ylabel('Average Path Length (hops)')
    plt.grid(True)
    plt.legend()
    plt.savefig(f"path_length_{config_key.replace(', ', '_')}.png")
    plt.show()


# --- Main script execution ---
if __name__ == "__main__":
    datasets = parse_results_file('experiment_results.txt')
    if datasets:
        for config, df in datasets.items():
            print(f"\nGenerating graphs for {config}...")
            create_plots(config, df)
        print("\nAll graphs have been generated and saved as PNG files.")