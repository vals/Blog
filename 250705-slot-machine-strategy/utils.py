import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List, Dict, Any

def encode_state(state: Tuple[str, ...], symbols: List[str]) -> int:
    index = 0
    num_symbols = len(symbols)
    num_reels = len(state)
    
    for i, symbol in enumerate(state):
        symbol_index = symbols.index(symbol)
        index += symbol_index * (num_symbols ** (num_reels - 1 - i))
    return index

def decode_state(index: int, symbols: List[str], num_reels: int) -> Tuple[str, ...]:
    state = []
    num_symbols = len(symbols)
    
    for i in range(num_reels):
        symbol_index = index // (num_symbols ** (num_reels - 1 - i))
        state.append(symbols[symbol_index])
        index = index % (num_symbols ** (num_reels - 1 - i))
    return tuple(state)

def print_q_table_summary(q_table: np.ndarray, top_n: int = 10):
    print(f"Q-table shape: {q_table.shape}")
    print(f"Max Q-value: {np.max(q_table):.4f}")
    print(f"Min Q-value: {np.min(q_table):.4f}")
    print(f"Mean Q-value: {np.mean(q_table):.4f}")
    
    max_q_states = np.unravel_index(np.argsort(q_table.ravel())[-top_n:], q_table.shape)
    print(f"\nTop {top_n} state-action pairs by Q-value:")
    for i in range(top_n):
        state_idx = max_q_states[0][-(i+1)]
        action_idx = max_q_states[1][-(i+1)]
        q_value = q_table[state_idx, action_idx]
        print(f"State {state_idx}, Action {action_idx}: Q = {q_value:.4f}")

def get_policy_from_q_table(q_table: np.ndarray) -> np.ndarray:
    return np.argmax(q_table, axis=1)

def plot_learning_curve(rewards: List[float], window_size: int = 100, filename: str = "learning_curve.png"):
    plt.figure(figsize=(12, 4))
    
    plt.subplot(1, 2, 1)
    plt.plot(rewards)
    plt.title('Episode Rewards')
    plt.xlabel('Episode')
    plt.ylabel('Reward')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    if len(rewards) >= window_size:
        moving_avg = np.convolve(rewards, np.ones(window_size)/window_size, mode='valid')
        plt.plot(moving_avg)
        plt.title(f'Moving Average Rewards (window={window_size})')
        plt.xlabel('Episode')
        plt.ylabel('Average Reward')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure to free memory
    print(f"Learning curve saved to {filename}")

def analyze_policy(policy: np.ndarray, game, actions: Dict[int, Dict[str, Any]]):
    print("Policy Analysis:")
    print("=" * 50)
    
    action_counts = {}
    for action in policy:
        action_counts[action] = action_counts.get(action, 0) + 1
    
    print("Action distribution:")
    for action, count in action_counts.items():
        percentage = (count / len(policy)) * 100
        description = actions[action]['description']
        cost = actions[action]['cost']
        print(f"Action {action} ({description}, cost {cost}): {count} states ({percentage:.1f}%)")
    
    print(f"\nMost common action: {max(action_counts, key=action_counts.get)}")
    print(f"Least common action: {min(action_counts, key=action_counts.get)}")

def plot_policy_distribution(policy: np.ndarray, actions: Dict[int, Dict[str, Any]], filename: str = "policy_distribution.png"):
    action_counts = {}
    for action in policy:
        action_counts[action] = action_counts.get(action, 0) + 1
    
    # Create bar chart
    plt.figure(figsize=(10, 6))
    action_labels = []
    counts = []
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FECA57']
    
    for action in sorted(action_counts.keys()):
        action_labels.append(f"Action {action}\n{actions[action]['description']}")
        counts.append(action_counts[action])
    
    bars = plt.bar(range(len(action_labels)), counts, color=colors[:len(action_labels)])
    plt.xlabel('Actions')
    plt.ylabel('Number of States')
    plt.title('Policy Action Distribution')
    plt.xticks(range(len(action_labels)), action_labels, rotation=45, ha='right')
    plt.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01,
                f'{count}\n({count/len(policy)*100:.1f}%)', 
                ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Policy distribution plot saved to {filename}")

def save_results(q_table: np.ndarray, rewards: List[float], policy: np.ndarray, filename: str = "results"):
    np.savez(f"{filename}.npz", 
             q_table=q_table, 
             rewards=np.array(rewards), 
             policy=policy)
    print(f"Results saved to {filename}.npz")

def load_results(filename: str = "results.npz"):
    data = np.load(filename)
    return data['q_table'], data['rewards'].tolist(), data['policy']