#!/usr/bin/env python3

"""
Simulate gameplay using the optimal expected value strategy.
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
from optimal_ev_qtable import OptimalEVQTable, OptimalEVAgent

class OptimalGameplaySimulator:
    def __init__(self, qtable_path: str = "optimal_qtable.npz"):
        self.qtable_path = qtable_path
        self.game = SlotMachine()
        self.qtable_builder = None
        self.optimal_agent = None
        self.balance_history = []
        self.episode_stats = []
        self.detailed_data = []  # For CSV export
        
    def load_optimal_qtable(self) -> bool:
        """Load or build the optimal Q-table."""
        try:
            if os.path.exists(self.qtable_path):
                # Load existing Q-table
                self.qtable_builder = OptimalEVQTable(self.game)
                self.qtable_builder.load_qtable(self.qtable_path)
                print(f"✅ Successfully loaded optimal Q-table from {self.qtable_path}")
            else:
                # Build Q-table if it doesn't exist
                print(f"❗ Q-table file not found: {self.qtable_path}")
                print("Building optimal Q-table (this may take a moment)...")
                self.qtable_builder = OptimalEVQTable(self.game)
                self.qtable_builder.build_optimal_qtable()
                self.qtable_builder.save_qtable(self.qtable_path)
                print(f"✅ Built and saved optimal Q-table to {self.qtable_path}")
            
            # Create optimal agent
            self.optimal_agent = OptimalEVAgent(self.game, self.qtable_builder)
            print(f"   Q-table shape: {self.qtable_builder.q_table.shape}")
            print(f"   Strategy: Mathematically optimal expected value decisions")
            return True
            
        except Exception as e:
            print(f"❌ Error with optimal Q-table: {e}")
            return False
    
    def simulate_episode(self, starting_balance: int, run_id: int = 0, episode_id: int = 0) -> Tuple[int, Dict]:
        """
        Simulate a single episode using the optimal strategy.
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
        
        # Play the episode using optimal policy
        max_steps = 20  # Same as other simulations
        for step in range(max_steps):
            # Record detailed data before taking action
            self._record_decision_data(
                run_id=run_id,
                episode_id=episode_id,
                balance=current_balance,
                step=step
            )
            
            # Choose action using optimal strategy
            action = self.optimal_agent.choose_action(self.game.current_state, self.game.respins_used)
            episode_actions.append(action)
            
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
    
    def _record_decision_data(self, run_id: int, episode_id: int, balance: int, step: int):
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
        
        # Get optimal action
        action = self.optimal_agent.choose_action(state, respins_used)
        
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
            print(f"🎰 Starting optimal strategy simulation:")
            print(f"   Initial balance: {initial_balance} Gold")
            print(f"   Maximum episodes: {max_episodes}")
            print(f"   Strategy: Mathematically optimal expected value decisions")
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
            print("🎯 OPTIMAL SIMULATION COMPLETED")
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
    
    def save_simulation_data(self, filename: str = "optimal_slot_machine_data.csv"):
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
        
        print(f"📊 Optimal simulation data saved to {filename}")
        print(f"   Total decisions recorded: {len(self.detailed_data)}")
    
    def simulate_multiple_sessions(self, num_runs: int = 25, initial_balance: int = 100, max_episodes: int = 500) -> List[List[int]]:
        """
        Run multiple independent simulation sessions.
        Returns list of balance histories for each run.
        """
        print(f"🎰 Running {num_runs} independent optimal strategy sessions:")
        print(f"   Initial balance: {initial_balance} Gold per session")
        print(f"   Maximum episodes: {max_episodes} per session")
        print(f"   Strategy: Mathematically optimal expected value decisions")
        print("-" * 50)
        
        all_balance_histories = []
        all_session_stats = []
        
        # Clear detailed data for fresh collection
        self.detailed_data = []
        
        # Run simulations with progress bar
        for run in tqdm(range(num_runs), desc="Running optimal simulations"):
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
        print("📊 AGGREGATE OPTIMAL STRATEGY RESULTS")
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
        
        # Strategy analysis
        if hasattr(self.optimal_agent, 'action_choices'):
            total_decisions = sum(self.optimal_agent.action_choices.values())
            if total_decisions > 0:
                cash_out_pct = (self.optimal_agent.action_choices['cash_out'] / total_decisions) * 100
                respin_pct = (self.optimal_agent.action_choices['respin'] / total_decisions) * 100
                print()
                print("DECISION ANALYSIS:")
                print(f"  Cash-out decisions: {cash_out_pct:.1f}%")
                print(f"  Respin decisions: {respin_pct:.1f}%")
                print(f"  Total decisions: {total_decisions}")
    
    def plot_multiple_balance_histories(self, all_balance_histories: List[List[int]], 
                                      all_session_stats: List[Dict], 
                                      filename: str = "optimal_balance_simulation.png"):
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
        
        ax1.set_title(f'Optimal Strategy Balance Progression - {len(all_balance_histories)} Independent Runs', 
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
        
        stats_text = f"""Optimal Strategy Statistics:
Total Runs: {len(all_balance_histories)}
Profitable: {profitable_runs} ({profitable_runs/len(all_balance_histories)*100:.1f}%)
Bankrupt: {bankruptcy_runs} ({bankruptcy_runs/len(all_balance_histories)*100:.1f}%)
Mean Final: {np.mean(final_balances):.1f} ± {np.std(final_balances):.1f}
Median Final: {np.median(final_balances):.1f}"""
        
        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"📊 Optimal strategy balance progressions saved to {filename}")

def main():
    """Main simulation function."""
    print("🎰 Optimal Slot Machine Strategy Simulation - Multiple Runs")
    print("=" * 60)
    print("Strategy: Mathematically optimal expected value decisions")
    print("=" * 60)
    
    # Initialize simulator
    simulator = OptimalGameplaySimulator()
    
    # Load optimal Q-table
    if not simulator.load_optimal_qtable():
        return
    
    print("\nRunning 25 independent 500-episode sessions...")
    print("This will show the performance of the mathematically optimal strategy.")
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
    with open("optimal_aggregate_results.txt", 'w') as f:
        f.write("OPTIMAL SLOT MACHINE STRATEGY - AGGREGATE SIMULATION RESULTS\n")
        f.write("=" * 60 + "\n\n")
        f.write("Strategy: Mathematically optimal expected value decisions\n\n")
        
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
    
    print("\n🎯 Optimal strategy simulation complete! Check the generated files:")
    print("   📊 optimal_balance_simulation.png - All 25 balance progression curves")
    print("   📄 optimal_aggregate_results.txt - Comprehensive statistics")
    print("   📊 optimal_slot_machine_data.csv - Detailed decision-by-decision data")
    print(f"\nThis analysis shows the performance of mathematically optimal play,")
    print(f"providing a theoretical benchmark for comparison with other strategies.")

if __name__ == "__main__":
    main()