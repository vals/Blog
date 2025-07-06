import json
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def load_empirical_data(filename):
    """Load empirical data from JSON export."""
    with open(filename, 'r') as f:
        data = json.load(f)
    
    # Extract balance history per run
    runs_data = defaultdict(list)
    
    # Group episodes by run
    for episode in data['episodes']:
        run_id = episode['runId']
        runs_data[run_id].append({
            'episodeId': episode['episodeId'],
            'endBalance': episode['endBalance'],
            'reward': episode['reward']
        })
    
    # Convert to balance curves
    empirical_runs = []
    for run_id, episodes in runs_data.items():
        # Sort episodes by episode ID
        episodes.sort(key=lambda x: x['episodeId'])
        
        balance_curve = []
        for episode in episodes:
            balance_curve.append(episode['endBalance'])
        
        empirical_runs.append(balance_curve)
    
    return empirical_runs, data['statistics']

def load_simulated_data(filename):
    """Load simulated Q-learning data from NPZ file."""
    try:
        data = np.load(filename)
        rewards = data['rewards']
        
        # Convert episode rewards to balance curves
        simulated_runs = []
        
        # Assume we want to create multiple runs from the reward data
        # Split the rewards into runs of reasonable length
        run_length = 50  # Episodes per run
        starting_balance = 100
        
        for start_idx in range(0, len(rewards), run_length):
            end_idx = min(start_idx + run_length, len(rewards))
            run_rewards = rewards[start_idx:end_idx]
            
            # Convert rewards to cumulative balance
            balance_curve = []
            current_balance = starting_balance
            
            for reward in run_rewards:
                current_balance += reward
                balance_curve.append(current_balance)
            
            simulated_runs.append(balance_curve)
            
            # Stop after reasonable number of runs for comparison
            if len(simulated_runs) >= 10:
                break
        
        return simulated_runs
    except Exception as e:
        print(f"Error loading simulated data: {e}")
        return []

def plot_comparison(empirical_runs, simulated_runs, save_filename='empirical_vs_simulated_comparison.png'):
    """Create comparison plot of empirical vs simulated runs."""
    plt.figure(figsize=(14, 8))
    
    # Plot empirical runs
    for i, run in enumerate(empirical_runs):
        episodes = range(1, len(run) + 1)
        plt.plot(episodes, run, color='blue', alpha=0.7, linewidth=2, 
                label='Human Player' if i == 0 else '')
    
    # Plot simulated runs
    for i, run in enumerate(simulated_runs):
        episodes = range(1, len(run) + 1)
        plt.plot(episodes, run, color='red', alpha=0.5, linewidth=1, 
                label='Q-Learning Agent' if i == 0 else '')
    
    plt.axhline(y=100, color='gray', linestyle='--', alpha=0.5, label='Starting Balance')
    
    plt.xlabel('Episode', fontsize=12)
    plt.ylabel('Balance (Gold)', fontsize=12)
    plt.title('Human Player vs Q-Learning Agent Performance Comparison', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    
    # Add summary statistics
    if empirical_runs:
        empirical_final_balances = [run[-1] for run in empirical_runs]
        avg_empirical = np.mean(empirical_final_balances)
        plt.text(0.02, 0.98, f'Human Avg Final Balance: {avg_empirical:.1f}', 
                transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='blue', alpha=0.1))
    
    if simulated_runs:
        simulated_final_balances = [run[-1] for run in simulated_runs]
        avg_simulated = np.mean(simulated_final_balances)
        plt.text(0.02, 0.90, f'Q-Learning Avg Final Balance: {avg_simulated:.1f}', 
                transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='red', alpha=0.1))
    
    plt.tight_layout()
    plt.savefig(save_filename, dpi=300, bbox_inches='tight')
    plt.show()
    print(f"Comparison plot saved to {save_filename}")

def analyze_performance(empirical_runs, simulated_runs):
    """Analyze and compare performance metrics."""
    print("PERFORMANCE ANALYSIS")
    print("=" * 50)
    
    if empirical_runs:
        empirical_final_balances = [run[-1] for run in empirical_runs]
        empirical_episode_counts = [len(run) for run in empirical_runs]
        
        print(f"HUMAN PLAYER:")
        print(f"  Number of runs: {len(empirical_runs)}")
        print(f"  Average episodes per run: {np.mean(empirical_episode_counts):.1f}")
        print(f"  Average final balance: {np.mean(empirical_final_balances):.1f}")
        print(f"  Best run final balance: {np.max(empirical_final_balances):.1f}")
        print(f"  Worst run final balance: {np.min(empirical_final_balances):.1f}")
        print(f"  Standard deviation: {np.std(empirical_final_balances):.1f}")
        
        # Calculate success rate (balance > 100)
        success_rate = len([b for b in empirical_final_balances if b > 100]) / len(empirical_final_balances)
        print(f"  Success rate (balance > 100): {success_rate:.1%}")
    
    if simulated_runs:
        simulated_final_balances = [run[-1] for run in simulated_runs]
        simulated_episode_counts = [len(run) for run in simulated_runs]
        
        print(f"\nQ-LEARNING AGENT:")
        print(f"  Number of runs: {len(simulated_runs)}")
        print(f"  Average episodes per run: {np.mean(simulated_episode_counts):.1f}")
        print(f"  Average final balance: {np.mean(simulated_final_balances):.1f}")
        print(f"  Best run final balance: {np.max(simulated_final_balances):.1f}")
        print(f"  Worst run final balance: {np.min(simulated_final_balances):.1f}")
        print(f"  Standard deviation: {np.std(simulated_final_balances):.1f}")
        
        # Calculate success rate (balance > 100)
        success_rate = len([b for b in simulated_final_balances if b > 100]) / len(simulated_final_balances)
        print(f"  Success rate (balance > 100): {success_rate:.1%}")
    
    if empirical_runs and simulated_runs:
        empirical_avg = np.mean([run[-1] for run in empirical_runs])
        simulated_avg = np.mean([run[-1] for run in simulated_runs])
        
        print(f"\nCOMPARISON:")
        if empirical_avg > simulated_avg:
            print(f"  Human player outperforms Q-learning by {empirical_avg - simulated_avg:.1f} gold on average")
        else:
            print(f"  Q-learning outperforms human player by {simulated_avg - empirical_avg:.1f} gold on average")

def main():
    # Load empirical data
    print("Loading empirical data...")
    try:
        empirical_runs, empirical_stats = load_empirical_data('slot_machine_data.json')
        print(f"Loaded {len(empirical_runs)} empirical runs")
    except Exception as e:
        print(f"Error loading empirical data: {e}")
        empirical_runs = []
    
    # Load simulated data
    print("Loading simulated data...")
    simulated_runs = []
    
    # Try different simulation result files
    sim_files = ['slot_machine_results.npz', 'enhanced_slot_machine_results.npz', 'advanced_slot_machine_results.npz']
    
    for sim_file in sim_files:
        try:
            runs = load_simulated_data(sim_file)
            if runs:
                simulated_runs.extend(runs)
                print(f"Loaded {len(runs)} simulated runs from {sim_file}")
                break
        except Exception as e:
            print(f"Could not load {sim_file}: {e}")
            continue
    
    if not simulated_runs:
        print("Warning: No simulated data loaded. Creating mock data for demonstration.")
        # Create some mock simulated runs for demonstration
        np.random.seed(42)
        for i in range(5):
            mock_run = []
            balance = 100
            for episode in range(30):
                # Mock Q-learning performance (slightly negative trend)
                reward = np.random.normal(-0.5, 3)
                balance += reward
                mock_run.append(balance)
            simulated_runs.append(mock_run)
    
    # Create comparison plot
    if empirical_runs or simulated_runs:
        plot_comparison(empirical_runs, simulated_runs)
        analyze_performance(empirical_runs, simulated_runs)
    else:
        print("No data available for comparison.")

if __name__ == "__main__":
    main()