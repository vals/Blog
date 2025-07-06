import json
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns

def load_and_analyze_empirical_data(filename):
    """Load empirical data and extract detailed metrics."""
    with open(filename, 'r') as f:
        data = json.load(f)
    
    # Extract detailed episode data
    episodes_by_run = defaultdict(list)
    all_episodes = []
    
    for episode in data['episodes']:
        run_id = episode['runId']
        episodes_by_run[run_id].append(episode)
        all_episodes.append(episode)
    
    # Extract decision patterns
    action_counts = defaultdict(int)
    for decision in data['decisions']:
        action_counts[decision['action']] += 1
    
    return episodes_by_run, all_episodes, action_counts, data

def create_comprehensive_comparison(empirical_data, simulated_data_file):
    """Create a comprehensive comparison plot."""
    episodes_by_run, all_episodes, action_counts, full_data = empirical_data
    
    # Load simulated data
    sim_data = np.load(simulated_data_file)
    sim_rewards = sim_data['rewards']
    
    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Balance progression comparison
    ax1 = plt.subplot(2, 3, (1, 2))
    
    # Plot empirical runs
    for run_id, episodes in episodes_by_run.items():
        episodes.sort(key=lambda x: x['episodeId'])
        balances = [ep['endBalance'] for ep in episodes]
        episode_nums = range(1, len(balances) + 1)
        plt.plot(episode_nums, balances, color='blue', alpha=0.6, linewidth=1.5)
    
    # Plot simulated data (convert to balance curves)
    starting_balance = 100
    run_length = 50
    for start_idx in range(0, min(len(sim_rewards), 500), run_length):
        end_idx = min(start_idx + run_length, len(sim_rewards))
        run_rewards = sim_rewards[start_idx:end_idx]
        
        balance_curve = []
        current_balance = starting_balance
        for reward in run_rewards:
            current_balance += reward
            balance_curve.append(current_balance)
        
        episode_nums = range(1, len(balance_curve) + 1)
        plt.plot(episode_nums, balance_curve, color='red', alpha=0.4, linewidth=1)
    
    plt.axhline(y=100, color='gray', linestyle='--', alpha=0.7, label='Starting Balance')
    plt.xlabel('Episode Number')
    plt.ylabel('Balance (Gold)')
    plt.title('Balance Progression: Human vs Q-Learning')
    plt.legend(['Human Runs', 'Q-Learning Runs', 'Starting Balance'])
    plt.grid(True, alpha=0.3)
    
    # 2. Final balance distribution
    ax2 = plt.subplot(2, 3, 3)
    
    empirical_final_balances = []
    for episodes in episodes_by_run.values():
        if episodes:
            episodes.sort(key=lambda x: x['episodeId'])
            empirical_final_balances.append(episodes[-1]['endBalance'])
    
    # Create simulated final balances
    simulated_final_balances = []
    for start_idx in range(0, min(len(sim_rewards), 500), run_length):
        end_idx = min(start_idx + run_length, len(sim_rewards))
        run_rewards = sim_rewards[start_idx:end_idx]
        final_balance = starting_balance + sum(run_rewards)
        simulated_final_balances.append(final_balance)
    
    plt.hist(empirical_final_balances, bins=15, alpha=0.7, color='blue', label='Human', density=True)
    plt.hist(simulated_final_balances, bins=15, alpha=0.7, color='red', label='Q-Learning', density=True)
    plt.xlabel('Final Balance')
    plt.ylabel('Density')
    plt.title('Final Balance Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 3. Action distribution (Human only)
    ax3 = plt.subplot(2, 3, 4)
    
    action_names = {
        0: 'Respin Blank',
        1: 'Respin Coin',
        2: 'Respin Coin Stack',
        3: 'Respin Snake',
        4: 'Respin Net',
        5: 'Respin x2',
        6: 'Respin Clover',
        7: 'Respin Crown',
        8: 'Cash Out'
    }
    
    actions = list(action_counts.keys())
    counts = list(action_counts.values())
    colors = ['skyblue' if action != 8 else 'lightcoral' for action in actions]
    
    plt.bar(range(len(actions)), counts, color=colors)
    plt.xlabel('Action')
    plt.ylabel('Frequency')
    plt.title('Human Action Distribution')
    plt.xticks(range(len(actions)), [action_names.get(a, f'Action {a}') for a in actions], rotation=45)
    plt.grid(True, alpha=0.3, axis='y')
    
    # 4. Episode rewards comparison
    ax4 = plt.subplot(2, 3, 5)
    
    empirical_rewards = [ep['reward'] for ep in all_episodes]
    
    plt.hist(empirical_rewards, bins=20, alpha=0.7, color='blue', label='Human', density=True)
    plt.hist(sim_rewards[:len(empirical_rewards)], bins=20, alpha=0.7, color='red', label='Q-Learning', density=True)
    plt.xlabel('Episode Reward')
    plt.ylabel('Density')
    plt.title('Episode Reward Distribution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 5. Performance metrics
    ax5 = plt.subplot(2, 3, 6)
    ax5.axis('off')
    
    # Calculate metrics
    human_avg_reward = np.mean(empirical_rewards)
    human_avg_final = np.mean(empirical_final_balances)
    human_success_rate = len([b for b in empirical_final_balances if b > 100]) / len(empirical_final_balances)
    
    sim_avg_reward = np.mean(sim_rewards[:len(empirical_rewards)])
    sim_avg_final = np.mean(simulated_final_balances)
    sim_success_rate = len([b for b in simulated_final_balances if b > 100]) / len(simulated_final_balances)
    
    metrics_text = f"""
PERFORMANCE COMPARISON

Human Player:
• Runs: {len(episodes_by_run)}
• Avg Episode Reward: {human_avg_reward:.2f}
• Avg Final Balance: {human_avg_final:.1f}
• Success Rate: {human_success_rate:.1%}
• Total Decisions: {len(full_data['decisions'])}

Q-Learning Agent:
• Runs: {len(simulated_final_balances)}
• Avg Episode Reward: {sim_avg_reward:.2f}
• Avg Final Balance: {sim_avg_final:.1f}
• Success Rate: {sim_success_rate:.1%}

Difference:
• Human vs AI Reward: {human_avg_reward - sim_avg_reward:.2f}
• Human vs AI Final: {human_avg_final - sim_avg_final:.1f}
    """
    
    ax5.text(0.1, 0.9, metrics_text, transform=ax5.transAxes, fontsize=10, 
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('detailed_empirical_vs_simulated.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Detailed comparison plot saved to 'detailed_empirical_vs_simulated.png'")

def main():
    # Load empirical data
    print("Loading empirical data...")
    empirical_data = load_and_analyze_empirical_data('slot_machine_data.json')
    
    # Create comprehensive comparison
    print("Creating detailed comparison...")
    create_comprehensive_comparison(empirical_data, 'slot_machine_results.npz')

if __name__ == "__main__":
    main()