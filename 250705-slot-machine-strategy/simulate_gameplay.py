#!/usr/bin/env python3

"""
Simulate gameplay using the learned Q-learning strategy.
Tracks balance progression over episodes and generates visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import os
import csv
import time
from tqdm import tqdm
from game_rules import SlotMachine
from q_learning import QLearningAgent

class GameplaySimulator:
    def __init__(self, model_path: str = "slot_machine_model.npz"):
        self.model_path = model_path
        self.game = SlotMachine()
        self.agent = None
        self.balance_history = []
        self.episode_stats = []
        self.detailed_data = []  # For CSV export
        
    def load_trained_model(self) -> bool:
        """Load the trained Q-learning model."""
        try:
            if not os.path.exists(self.model_path):
                print(f"❌ Model file not found: {self.model_path}")
                print("Please run main.py first to train the model.")
                return False
            
            # Initialize agent with same parameters as training
            self.agent = QLearningAgent(
                state_space_size=self.game.get_abstract_state_space_size(),
                action_space_size=self.game.get_action_space_size(),
                learning_rate=0.25,
                discount_factor=0.99,
                epsilon=0.0,  # No exploration during simulation
                epsilon_decay=1.0,  # No decay needed
                epsilon_min=0.0
            )
            
            # Load the trained model
            self.agent.load_model(self.model_path)
            print(f"✅ Successfully loaded trained model from {self.model_path}")
            print(f"   Q-table shape: {self.agent.q_table.shape}")
            print(f"   Non-zero Q-values: {np.count_nonzero(self.agent.q_table)}")
            return True
            
        except Exception as e:
            print(f"❌ Error loading model: {e}")
            return False
    
    def simulate_episode(self, starting_balance: int, run_id: int = 0, episode_id: int = 0) -> Tuple[int, Dict]:
        """
        Simulate a single episode using the learned policy.
        Returns final balance and episode statistics.
        """
        if starting_balance < 1:
            return starting_balance, {
                'initial_balance': starting_balance,
                'final_balance': starting_balance,
                'total_cost': 0,
                'payout': 0,
                'net_result': 0,
                'actions_taken': 0,
                'episode_length': 0,
                'insufficient_funds': True
            }
        
        # Reset game for new episode
        state = self.game.reset()
        initial_balance = starting_balance
        episode_actions = []
        current_balance = starting_balance
        
        # Play the episode using learned policy
        max_steps = 20  # Same as training
        for step in range(max_steps):
            # Get current state and valid actions
            state_index = self.game.get_current_abstract_state_index()
            valid_actions = self.game.get_valid_actions()
            
            # Choose action using learned policy (no exploration)
            action = self.agent.choose_action(state_index, valid_actions)
            episode_actions.append(action)
            
            # Record detailed data before taking action
            self._record_decision_data(
                run_id=run_id,
                episode_id=episode_id,
                action=action,
                balance=current_balance,
                step=step
            )
            
            # Take action
            reward, next_state, done = self.game.step(action)
            
            # Update balance after action
            current_balance = starting_balance + reward
            
            if done:
                break
        
        # Calculate episode results
        total_cost = self.game.total_cost
        payout = total_cost + reward  # reward is net (payout - cost)
        net_result = reward
        final_balance = starting_balance + net_result
        
        episode_stats = {
            'initial_balance': initial_balance,
            'final_balance': final_balance,
            'total_cost': total_cost,
            'payout': payout,
            'net_result': net_result,
            'actions_taken': len(episode_actions),
            'episode_length': step + 1,
            'insufficient_funds': False,
            'actions': episode_actions
        }
        
        return final_balance, episode_stats
    
    def _record_decision_data(self, run_id: int, episode_id: int, action: int, balance: int, step: int):
        """Record detailed decision data for CSV export matching empirical data format."""
        # Get current game state
        state = self.game.current_state
        respins_used = self.game.respins_used
        total_cost = self.game.total_cost
        
        # Get current payout (before action)
        current_payout = self.game._calculate_payout(state, 1)
        
        # Get symbol counts
        symbol_counts = {symbol: 0 for symbol in self.game.symbols}
        for symbol in state:
            symbol_counts[symbol] += 1
        
        # Get valid actions
        valid_actions = self.game.get_valid_actions()
        valid_actions_str = ";".join(map(str, valid_actions))
        
        # Create decision record
        decision_id = len(self.detailed_data) + 1
        timestamp = int(time.time() * 1000)  # milliseconds
        
        record = {
            'decisionId': decision_id,
            'runId': run_id,
            'episodeId': episode_id,
            'timestamp': timestamp,
            'action': action,
            'actionDescription': self.game.actions[action]['description'],
            'respinsUsed': respins_used,
            'balance': balance,
            'currentPayout': current_payout,
            'totalCost': total_cost,
            'blank': symbol_counts['Blank'],
            'coin': symbol_counts['Coin'],
            'coinStack': symbol_counts['Coin Stack'],
            'snake': symbol_counts['Snake'],
            'net': symbol_counts['Net'],
            'x2': symbol_counts['x2'],
            'clover': symbol_counts['Clover'],
            'crown': symbol_counts['Crown'],
            'reel0': state[0],
            'reel1': state[1],
            'reel2': state[2],
            'reel3': state[3],
            'validActions': valid_actions_str
        }
        
        self.detailed_data.append(record)
    
    def simulate_session(self, initial_balance: int = 100, max_episodes: int = 500, verbose: bool = True, run_id: int = 0) -> Dict:
        """
        Simulate a full gameplay session until balance reaches 0 or max episodes.
        """
        if verbose:
            print(f"🎰 Starting gameplay simulation:")
            print(f"   Initial balance: {initial_balance} Gold")
            print(f"   Maximum episodes: {max_episodes}")
            print(f"   Stop condition: Balance ≤ 0 or {max_episodes} episodes")
            print("-" * 50)
        
        balance = initial_balance
        self.balance_history = [balance]
        self.episode_stats = []
        
        episode = 0
        while episode < max_episodes and balance > 0:
            # Simulate episode
            new_balance, stats = self.simulate_episode(balance, run_id, episode + 1)
            
            # Update tracking
            balance = new_balance
            self.balance_history.append(balance)
            self.episode_stats.append(stats)
            episode += 1
            
            # Progress reporting every 50 episodes
            if verbose and (episode % 50 == 0 or balance <= 0):
                print(f"Episode {episode:3d}: Balance = {balance:3d} Gold "
                      f"(Net: {stats['net_result']:+3d})")
        
        # Session summary
        session_stats = self._calculate_session_stats(initial_balance, episode)
        
        if verbose:
            print("\n" + "="*50)
            print("🎯 SIMULATION COMPLETED")
            print("="*50)
            print(f"Episodes played: {episode}")
            print(f"Final balance: {balance} Gold")
            print(f"Net result: {balance - initial_balance:+d} Gold")
            print(f"Session outcome: {'💰 Profit' if balance > initial_balance else '💸 Loss' if balance < initial_balance else '🤝 Break-even'}")
        
        return session_stats
    
    def _calculate_session_stats(self, initial_balance: int, episodes_played: int) -> Dict:
        """Calculate comprehensive session statistics."""
        if not self.episode_stats:
            return {}
        
        # Basic stats
        final_balance = self.balance_history[-1]
        net_result = final_balance - initial_balance
        
        # Episode results
        profits = [ep['net_result'] for ep in self.episode_stats if ep['net_result'] > 0]
        losses = [ep['net_result'] for ep in self.episode_stats if ep['net_result'] < 0]
        breakevens = [ep for ep in self.episode_stats if ep['net_result'] == 0]
        
        # Streaks
        max_balance = max(self.balance_history)
        min_balance = min(self.balance_history)
        
        # Performance metrics
        win_rate = len(profits) / len(self.episode_stats) * 100 if self.episode_stats else 0
        avg_profit = np.mean(profits) if profits else 0
        avg_loss = np.mean(losses) if losses else 0
        avg_episode_result = np.mean([ep['net_result'] for ep in self.episode_stats])
        
        return {
            'initial_balance': initial_balance,
            'final_balance': final_balance,
            'net_result': net_result,
            'episodes_played': episodes_played,
            'max_balance_reached': max_balance,
            'min_balance_reached': min_balance,
            'winning_episodes': len(profits),
            'losing_episodes': len(losses),
            'breakeven_episodes': len(breakevens),
            'win_rate': win_rate,
            'avg_profit_per_win': avg_profit,
            'avg_loss_per_loss': avg_loss,
            'avg_episode_result': avg_episode_result,
            'total_wagered': sum(ep['total_cost'] for ep in self.episode_stats),
            'total_payout': sum(ep['payout'] for ep in self.episode_stats)
        }
    
    def plot_balance_history(self, session_stats: Dict, filename: str = "balance_simulation.png"):
        """Create visualization of balance progression."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
        
        # Main balance progression plot
        episodes = range(len(self.balance_history))
        ax1.plot(episodes, self.balance_history, linewidth=2, color='#2E86C1', alpha=0.8)
        ax1.fill_between(episodes, self.balance_history, alpha=0.3, color='#2E86C1')
        
        # Add horizontal lines for initial balance and zero
        ax1.axhline(y=session_stats['initial_balance'], color='green', linestyle='--', alpha=0.7, label='Initial Balance')
        ax1.axhline(y=0, color='red', linestyle='--', alpha=0.7, label='Bankruptcy')
        
        # Highlight max and min points
        max_ep = self.balance_history.index(session_stats['max_balance_reached'])
        min_ep = self.balance_history.index(session_stats['min_balance_reached'])
        ax1.scatter(max_ep, session_stats['max_balance_reached'], color='green', s=100, zorder=5, label='Peak Balance')
        ax1.scatter(min_ep, session_stats['min_balance_reached'], color='red', s=100, zorder=5, label='Lowest Balance')
        
        ax1.set_title('Balance Progression Over Episodes', fontsize=14, fontweight='bold')
        ax1.set_xlabel('Episode')
        ax1.set_ylabel('Balance (Gold)')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Episode results histogram
        episode_results = [ep['net_result'] for ep in self.episode_stats]
        ax2.hist(episode_results, bins=30, alpha=0.7, color='#28B463', edgecolor='black')
        ax2.axvline(x=0, color='red', linestyle='--', alpha=0.7)
        ax2.set_title('Distribution of Episode Results', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Net Result per Episode (Gold)')
        ax2.set_ylabel('Frequency')
        ax2.grid(True, alpha=0.3)
        
        # Add statistics text box
        stats_text = f"""Session Statistics:
Episodes: {session_stats['episodes_played']}
Final Balance: {session_stats['final_balance']} Gold
Net Result: {session_stats['net_result']:+d} Gold
Win Rate: {session_stats['win_rate']:.1f}%
Peak Balance: {session_stats['max_balance_reached']} Gold
Avg Episode: {session_stats['avg_episode_result']:+.2f} Gold"""
        
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Balance progression saved to {filename}")
    
    def save_detailed_stats(self, session_stats: Dict, filename: str = "simulation_results.txt"):
        """Save detailed session statistics to text file."""
        with open(filename, 'w') as f:
            f.write("SLOT MACHINE STRATEGY SIMULATION RESULTS\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("SESSION OVERVIEW:\n")
            f.write(f"Initial Balance: {session_stats['initial_balance']} Gold\n")
            f.write(f"Final Balance: {session_stats['final_balance']} Gold\n")
            f.write(f"Net Result: {session_stats['net_result']:+d} Gold\n")
            f.write(f"Episodes Played: {session_stats['episodes_played']}\n")
            f.write(f"Peak Balance: {session_stats['max_balance_reached']} Gold\n")
            f.write(f"Lowest Balance: {session_stats['min_balance_reached']} Gold\n\n")
            
            f.write("PERFORMANCE METRICS:\n")
            f.write(f"Win Rate: {session_stats['win_rate']:.2f}%\n")
            f.write(f"Winning Episodes: {session_stats['winning_episodes']}\n")
            f.write(f"Losing Episodes: {session_stats['losing_episodes']}\n")
            f.write(f"Break-even Episodes: {session_stats['breakeven_episodes']}\n")
            f.write(f"Average Episode Result: {session_stats['avg_episode_result']:+.3f} Gold\n")
            f.write(f"Average Profit per Win: {session_stats['avg_profit_per_win']:.2f} Gold\n")
            f.write(f"Average Loss per Loss: {session_stats['avg_loss_per_loss']:.2f} Gold\n\n")
            
            f.write("WAGERING SUMMARY:\n")
            f.write(f"Total Wagered: {session_stats['total_wagered']} Gold\n")
            f.write(f"Total Payout: {session_stats['total_payout']} Gold\n")
            f.write(f"Return Rate: {session_stats['total_payout']/session_stats['total_wagered']*100:.2f}%\n")
        
        print(f"📄 Detailed statistics saved to {filename}")
    
    def save_simulation_data(self, filename: str = "simulated_slot_machine_data.csv"):
        """Save detailed simulation data to CSV in format matching empirical data."""
        if not self.detailed_data:
            print("⚠️  No detailed data to save. Run simulations first.")
            return
        
        # Define CSV headers matching empirical data format
        headers = [
            'decisionId', 'runId', 'episodeId', 'timestamp', 'action', 'actionDescription',
            'respinsUsed', 'balance', 'currentPayout', 'totalCost', 'blank', 'coin',
            'coinStack', 'snake', 'net', 'x2', 'clover', 'crown', 'reel0', 'reel1',
            'reel2', 'reel3', 'validActions'
        ]
        
        with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            
            for row in self.detailed_data:
                writer.writerow(row)
        
        print(f"📊 Simulation data saved to {filename}")
        print(f"   Total decisions recorded: {len(self.detailed_data)}")
    
    def simulate_multiple_sessions(self, num_runs: int = 25, initial_balance: int = 100, max_episodes: int = 500) -> List[List[int]]:
        """
        Run multiple independent simulation sessions.
        Returns list of balance histories for each run.
        """
        print(f"🎰 Running {num_runs} independent simulation sessions:")
        print(f"   Initial balance: {initial_balance} Gold per session")
        print(f"   Maximum episodes: {max_episodes} per session")
        print("-" * 50)
        
        all_balance_histories = []
        all_session_stats = []
        
        # Clear detailed data for fresh collection
        self.detailed_data = []
        
        # Run simulations with progress bar
        for run in tqdm(range(num_runs), desc="Running simulations"):
            # Reset for new run
            self.balance_history = []
            self.episode_stats = []
            
            # Run single session
            session_stats = self.simulate_session(
                initial_balance=initial_balance,
                max_episodes=max_episodes,
                verbose=False,  # Suppress individual session output
                run_id=run + 1  # 1-indexed run IDs
            )
            
            # Store results
            all_balance_histories.append(self.balance_history.copy())
            all_session_stats.append(session_stats)
            
            # Brief progress update
            if (run + 1) % 5 == 0:
                recent_finals = [stats['final_balance'] for stats in all_session_stats[-5:]]
                avg_final = np.mean(recent_finals)
                tqdm.write(f"Runs {run-4}-{run+1}: Avg final balance = {avg_final:.1f} Gold")
        
        # Calculate aggregate statistics
        self._print_aggregate_stats(all_session_stats, initial_balance)
        
        return all_balance_histories, all_session_stats
    
    def _print_aggregate_stats(self, all_session_stats: List[Dict], initial_balance: int):
        """Print aggregate statistics across all runs."""
        print("\n" + "="*60)
        print("📊 AGGREGATE RESULTS ACROSS ALL RUNS")
        print("="*60)
        
        final_balances = [stats['final_balance'] for stats in all_session_stats]
        net_results = [stats['net_result'] for stats in all_session_stats]
        episodes_played = [stats['episodes_played'] for stats in all_session_stats]
        
        # Basic statistics
        print(f"Total runs: {len(all_session_stats)}")
        print(f"Initial balance per run: {initial_balance} Gold")
        print()
        
        # Final balance statistics
        print("FINAL BALANCE STATISTICS:")
        print(f"  Mean final balance: {np.mean(final_balances):.2f} ± {np.std(final_balances):.2f} Gold")
        print(f"  Median final balance: {np.median(final_balances):.2f} Gold")
        print(f"  Min final balance: {np.min(final_balances)} Gold")
        print(f"  Max final balance: {np.max(final_balances)} Gold")
        print()
        
        # Performance statistics
        profitable_runs = sum(1 for balance in final_balances if balance > initial_balance)
        bankruptcy_runs = sum(1 for balance in final_balances if balance == 0)
        
        print("PERFORMANCE STATISTICS:")
        print(f"  Profitable runs: {profitable_runs}/{len(all_session_stats)} ({profitable_runs/len(all_session_stats)*100:.1f}%)")
        print(f"  Bankrupt runs: {bankruptcy_runs}/{len(all_session_stats)} ({bankruptcy_runs/len(all_session_stats)*100:.1f}%)")
        print(f"  Break-even runs: {len(all_session_stats) - profitable_runs - bankruptcy_runs}")
        print()
        
        # Episode statistics
        print("EPISODE STATISTICS:")
        print(f"  Mean episodes played: {np.mean(episodes_played):.1f} ± {np.std(episodes_played):.1f}")
        print(f"  Median episodes played: {np.median(episodes_played):.1f}")
        print()
        
        # Net result statistics
        print("NET RESULT STATISTICS:")
        print(f"  Mean net result: {np.mean(net_results):+.2f} ± {np.std(net_results):.2f} Gold")
        print(f"  Total net result: {np.sum(net_results):+.0f} Gold")
        print(f"  Best single run: {np.max(net_results):+.0f} Gold")
        print(f"  Worst single run: {np.min(net_results):+.0f} Gold")
    
    def plot_multiple_balance_histories(self, all_balance_histories: List[List[int]], 
                                      all_session_stats: List[Dict], 
                                      filename: str = "multiple_balance_simulation.png"):
        """Create visualization showing all balance progressions."""
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 12))
        
        # Plot all individual balance curves
        final_balances = []
        max_episodes = max(len(history) for history in all_balance_histories)
        
        for i, balance_history in enumerate(all_balance_histories):
            episodes = range(len(balance_history))
            final_balance = balance_history[-1]
            final_balances.append(final_balance)
            
            # Color coding: green for profit, red for loss, blue for break-even
            initial_balance = all_session_stats[0]['initial_balance']
            if final_balance > initial_balance:
                color = 'green'
                alpha = 0.6
            elif final_balance == 0:
                color = 'red'
                alpha = 0.8
            else:
                color = 'gray'
                alpha = 0.4
            
            ax1.plot(episodes, balance_history, color=color, alpha=alpha, linewidth=1)
        
        # Calculate and plot average curve
        # Pad all histories to same length for averaging
        padded_histories = []
        for history in all_balance_histories:
            padded = history + [history[-1]] * (max_episodes - len(history))
            padded_histories.append(padded)
        
        avg_balance = np.mean(padded_histories, axis=0)
        episodes_avg = range(len(avg_balance))
        ax1.plot(episodes_avg, avg_balance, color='black', linewidth=3, label='Average Balance')
        
        # Add reference lines
        initial_balance = all_session_stats[0]['initial_balance']
        ax1.axhline(y=initial_balance, color='blue', linestyle='--', alpha=0.7, label='Initial Balance')
        ax1.axhline(y=0, color='red', linestyle='--', alpha=0.7, label='Bankruptcy')
        
        ax1.set_title(f'Balance Progression - {len(all_balance_histories)} Independent Runs', 
                     fontsize=14, fontweight='bold')
        ax1.set_xlabel('Episode')
        ax1.set_ylabel('Balance (Gold)')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # Final balance distribution
        ax2.hist(final_balances, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax2.axvline(x=initial_balance, color='blue', linestyle='--', alpha=0.7, label='Initial Balance')
        ax2.axvline(x=np.mean(final_balances), color='red', linestyle='-', linewidth=2, label='Mean Final Balance')
        ax2.set_title('Distribution of Final Balances', fontsize=14, fontweight='bold')
        ax2.set_xlabel('Final Balance (Gold)')
        ax2.set_ylabel('Frequency')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # Add aggregate statistics text box
        profitable_runs = sum(1 for balance in final_balances if balance > initial_balance)
        bankruptcy_runs = sum(1 for balance in final_balances if balance == 0)
        
        stats_text = f"""Aggregate Statistics:
Total Runs: {len(all_balance_histories)}
Profitable: {profitable_runs} ({profitable_runs/len(all_balance_histories)*100:.1f}%)
Bankrupt: {bankruptcy_runs} ({bankruptcy_runs/len(all_balance_histories)*100:.1f}%)
Mean Final: {np.mean(final_balances):.1f} ± {np.std(final_balances):.1f}
Median Final: {np.median(final_balances):.1f}"""
        
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Multiple balance progressions saved to {filename}")

def main():
    """Main simulation function."""
    print("🎰 Slot Machine Strategy Simulation - Multiple Runs")
    print("=" * 60)
    
    # Initialize simulator
    simulator = GameplaySimulator()
    
    # Load trained model
    if not simulator.load_trained_model():
        return
    
    print("\nRunning 25 independent 500-episode sessions...")
    print("This will show the variability and robustness of the learned strategy.")
    print()
    
    # Run multiple simulations
    all_balance_histories, all_session_stats = simulator.simulate_multiple_sessions(
        num_runs=25,
        initial_balance=100,
        max_episodes=500
    )
    
    # Generate comprehensive visualization
    simulator.plot_multiple_balance_histories(all_balance_histories, all_session_stats)
    
    # Save detailed simulation data to CSV
    simulator.save_simulation_data()
    
    # Also save aggregate statistics
    with open("aggregate_simulation_results.txt", 'w') as f:
        f.write("SLOT MACHINE STRATEGY - AGGREGATE SIMULATION RESULTS\n")
        f.write("=" * 60 + "\n\n")
        
        final_balances = [stats['final_balance'] for stats in all_session_stats]
        net_results = [stats['net_result'] for stats in all_session_stats]
        episodes_played = [stats['episodes_played'] for stats in all_session_stats]
        
        f.write(f"Total independent runs: {len(all_session_stats)}\n")
        f.write(f"Episodes per run: 500 (or until bankruptcy)\n")
        f.write(f"Initial balance per run: 100 Gold\n\n")
        
        f.write("FINAL BALANCE STATISTICS:\n")
        f.write(f"Mean: {np.mean(final_balances):.2f} ± {np.std(final_balances):.2f} Gold\n")
        f.write(f"Median: {np.median(final_balances):.2f} Gold\n")
        f.write(f"Range: {np.min(final_balances)} - {np.max(final_balances)} Gold\n\n")
        
        profitable_runs = sum(1 for balance in final_balances if balance > 100)
        bankruptcy_runs = sum(1 for balance in final_balances if balance == 0)
        
        f.write("PERFORMANCE SUMMARY:\n")
        f.write(f"Profitable runs: {profitable_runs}/25 ({profitable_runs/25*100:.1f}%)\n")
        f.write(f"Bankrupt runs: {bankruptcy_runs}/25 ({bankruptcy_runs/25*100:.1f}%)\n")
        f.write(f"Break-even runs: {25 - profitable_runs - bankruptcy_runs}/25\n\n")
        
        f.write("NET RESULT STATISTICS:\n")
        f.write(f"Mean per run: {np.mean(net_results):+.2f} ± {np.std(net_results):.2f} Gold\n")
        f.write(f"Total across all runs: {np.sum(net_results):+.0f} Gold\n")
        f.write(f"Best single run: {np.max(net_results):+.0f} Gold\n")
        f.write(f"Worst single run: {np.min(net_results):+.0f} Gold\n")
    
    print("\n🎯 Multi-run simulation complete! Check the generated files:")
    print("   📊 multiple_balance_simulation.png - All 25 balance progression curves")
    print("   📄 aggregate_simulation_results.txt - Comprehensive statistics")
    print("   📊 simulated_slot_machine_data.csv - Detailed decision-by-decision data")
    print(f"\nThis analysis shows the robustness of the learned strategy across")
    print(f"multiple independent sessions, revealing both the typical performance")
    print(f"and the variability you can expect from the Q-learning policy.")

if __name__ == "__main__":
    main()