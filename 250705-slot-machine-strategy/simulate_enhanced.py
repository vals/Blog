#!/usr/bin/env python3

"""
Simulate gameplay using the enhanced Q-learning strategy.
Compare enhanced vs baseline performance in realistic gameplay.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import os
from tqdm import tqdm
from game_rules import SlotMachine
from enhanced_q_learning import EnhancedQLearningAgent
from q_learning import QLearningAgent

class EnhancedGameplaySimulator:
    def __init__(self, enhanced_model_path: str = "enhanced_slot_machine_model.npz",
                 baseline_model_path: str = "slot_machine_model.npz"):
        self.enhanced_model_path = enhanced_model_path
        self.baseline_model_path = baseline_model_path
        self.game = SlotMachine()
        self.enhanced_agent = None
        self.baseline_agent = None
        
    def load_models(self) -> bool:
        """Load both enhanced and baseline models."""
        success = True
        
        # Load enhanced model
        try:
            if os.path.exists(self.enhanced_model_path):
                self.enhanced_agent = EnhancedQLearningAgent(
                    state_space_size=self.game.get_abstract_state_space_size(),
                    action_space_size=self.game.get_action_space_size(),
                    epsilon=0.0  # No exploration during simulation
                )
                self.enhanced_agent.load_model(self.enhanced_model_path)
                print(f"✅ Enhanced model loaded from {self.enhanced_model_path}")
            else:
                print(f"❌ Enhanced model not found: {self.enhanced_model_path}")
                success = False
        except Exception as e:
            print(f"❌ Error loading enhanced model: {e}")
            success = False
        
        # Load baseline model
        try:
            if os.path.exists(self.baseline_model_path):
                self.baseline_agent = QLearningAgent(
                    state_space_size=self.game.get_abstract_state_space_size(),
                    action_space_size=self.game.get_action_space_size(),
                    epsilon=0.0  # No exploration during simulation
                )
                self.baseline_agent.load_model(self.baseline_model_path)
                print(f"✅ Baseline model loaded from {self.baseline_model_path}")
            else:
                print(f"❌ Baseline model not found: {self.baseline_model_path}")
                success = False
        except Exception as e:
            print(f"❌ Error loading baseline model: {e}")
            success = False
        
        return success
    
    def simulate_episode(self, agent, starting_balance: int) -> Tuple[int, Dict]:
        """Simulate a single episode using the given agent."""
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
        
        # Play the episode using the agent's policy
        max_steps = 20
        for step in range(max_steps):
            state_index = self.game.get_current_abstract_state_index()
            valid_actions = self.game.get_valid_actions()
            
            # Choose action using the agent's policy
            action = agent.choose_action(state_index, valid_actions)
            episode_actions.append(action)
            
            # Take action
            reward, next_state, done = self.game.step(action)
            
            if done:
                break
        
        # Calculate episode results
        total_cost = self.game.total_cost
        payout = total_cost + reward
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
    
    def simulate_agent_session(self, agent, agent_name: str, initial_balance: int = 100, 
                              max_episodes: int = 500, verbose: bool = False) -> Tuple[List[int], List[Dict]]:
        """Simulate a full session for one agent."""
        if verbose:
            print(f"🎰 Simulating {agent_name} session:")
            print(f"   Initial balance: {initial_balance} Gold")
            print(f"   Maximum episodes: {max_episodes}")
        
        balance = initial_balance
        balance_history = [balance]
        episode_stats = []
        
        episode = 0
        while episode < max_episodes and balance > 0:
            # Simulate episode
            new_balance, stats = self.simulate_episode(agent, balance)
            
            # Update tracking
            balance = new_balance
            balance_history.append(balance)
            episode_stats.append(stats)
            episode += 1
            
            # Progress reporting
            if verbose and (episode % 100 == 0 or balance <= 0):
                print(f"  Episode {episode:3d}: Balance = {balance:3d} Gold")
        
        if verbose:
            print(f"  Final result: {balance} Gold after {episode} episodes")
        
        return balance_history, episode_stats
    
    def compare_agents(self, num_runs: int = 25, initial_balance: int = 100, 
                      max_episodes: int = 500) -> Dict:
        """Compare enhanced vs baseline agents across multiple runs."""
        print(f"🎯 Comparing Enhanced vs Baseline Agents")
        print(f"   Running {num_runs} independent sessions per agent")
        print(f"   Initial balance: {initial_balance} Gold per session")
        print(f"   Maximum episodes: {max_episodes} per session")
        print("-" * 60)
        
        enhanced_results = []
        baseline_results = []
        
        # Run simulations for both agents
        for run in tqdm(range(num_runs), desc="Running comparisons"):
            # Enhanced agent
            enhanced_history, enhanced_stats = self.simulate_agent_session(
                self.enhanced_agent, "Enhanced", initial_balance, max_episodes, verbose=False
            )
            enhanced_results.append({
                'balance_history': enhanced_history,
                'episode_stats': enhanced_stats,
                'final_balance': enhanced_history[-1],
                'episodes_played': len(enhanced_stats)
            })
            
            # Baseline agent
            baseline_history, baseline_stats = self.simulate_agent_session(
                self.baseline_agent, "Baseline", initial_balance, max_episodes, verbose=False
            )
            baseline_results.append({
                'balance_history': baseline_history,
                'episode_stats': baseline_stats,
                'final_balance': baseline_history[-1],
                'episodes_played': len(baseline_stats)
            })
        
        # Analyze results
        comparison_stats = self._analyze_comparison(enhanced_results, baseline_results, initial_balance)
        
        # Create comparison visualization
        self._plot_agent_comparison(enhanced_results, baseline_results, comparison_stats)
        
        return comparison_stats
    
    def _analyze_comparison(self, enhanced_results: List[Dict], baseline_results: List[Dict], 
                           initial_balance: int) -> Dict:
        """Analyze and compare the results from both agents."""
        enhanced_finals = [r['final_balance'] for r in enhanced_results]
        baseline_finals = [r['final_balance'] for r in baseline_results]
        
        enhanced_episodes = [r['episodes_played'] for r in enhanced_results]
        baseline_episodes = [r['episodes_played'] for r in baseline_results]
        
        # Performance statistics
        enhanced_mean = np.mean(enhanced_finals)
        enhanced_std = np.std(enhanced_finals)
        baseline_mean = np.mean(baseline_finals)
        baseline_std = np.std(baseline_finals)
        
        # Win/loss statistics
        enhanced_profitable = sum(1 for f in enhanced_finals if f > initial_balance)
        enhanced_bankrupt = sum(1 for f in enhanced_finals if f == 0)
        baseline_profitable = sum(1 for f in baseline_finals if f > initial_balance)
        baseline_bankrupt = sum(1 for f in baseline_finals if f == 0)
        
        # Calculate improvement
        improvement = (enhanced_mean - baseline_mean) / abs(baseline_mean) * 100 if baseline_mean != 0 else 0
        
        comparison_stats = {
            'enhanced_mean': enhanced_mean,
            'enhanced_std': enhanced_std,
            'enhanced_profitable': enhanced_profitable,
            'enhanced_bankrupt': enhanced_bankrupt,
            'baseline_mean': baseline_mean,
            'baseline_std': baseline_std,
            'baseline_profitable': baseline_profitable,
            'baseline_bankrupt': baseline_bankrupt,
            'improvement_percent': improvement,
            'num_runs': len(enhanced_results)
        }
        
        # Print comparison
        print("\n" + "="*60)
        print("📊 AGENT COMPARISON RESULTS")
        print("="*60)
        
        print(f"Enhanced Agent:")
        print(f"  Mean final balance: {enhanced_mean:.2f} ± {enhanced_std:.2f} Gold")
        print(f"  Profitable runs: {enhanced_profitable}/{len(enhanced_results)} ({enhanced_profitable/len(enhanced_results)*100:.1f}%)")
        print(f"  Bankrupt runs: {enhanced_bankrupt}/{len(enhanced_results)} ({enhanced_bankrupt/len(enhanced_results)*100:.1f}%)")
        print(f"  Mean episodes: {np.mean(enhanced_episodes):.1f}")
        
        print(f"\nBaseline Agent:")
        print(f"  Mean final balance: {baseline_mean:.2f} ± {baseline_std:.2f} Gold")
        print(f"  Profitable runs: {baseline_profitable}/{len(baseline_results)} ({baseline_profitable/len(baseline_results)*100:.1f}%)")
        print(f"  Bankrupt runs: {baseline_bankrupt}/{len(baseline_results)} ({baseline_bankrupt/len(baseline_results)*100:.1f}%)")
        print(f"  Mean episodes: {np.mean(baseline_episodes):.1f}")
        
        print(f"\n🎯 IMPROVEMENT: {improvement:+.2f}%")
        
        if improvement > 5:
            print("🚀 Significant improvement with enhanced agent!")
        elif improvement > 0:
            print("✅ Modest improvement with enhanced agent.")
        else:
            print("⚠️  Enhanced agent performed worse than baseline.")
        
        return comparison_stats
    
    def _plot_agent_comparison(self, enhanced_results: List[Dict], baseline_results: List[Dict], 
                              comparison_stats: Dict):
        """Create comprehensive visualization comparing both agents."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        
        # Balance progression comparison
        max_episodes = max(
            max(len(r['balance_history']) for r in enhanced_results),
            max(len(r['balance_history']) for r in baseline_results)
        )
        
        # Plot individual curves
        for result in enhanced_results:
            episodes = range(len(result['balance_history']))
            ax1.plot(episodes, result['balance_history'], color='red', alpha=0.3, linewidth=1)
        
        for result in baseline_results:
            episodes = range(len(result['balance_history']))
            ax1.plot(episodes, result['balance_history'], color='blue', alpha=0.3, linewidth=1)
        
        # Plot average curves
        enhanced_avg = self._calculate_average_curve(enhanced_results, max_episodes)
        baseline_avg = self._calculate_average_curve(baseline_results, max_episodes)
        
        ax1.plot(range(len(enhanced_avg)), enhanced_avg, color='red', linewidth=3, label='Enhanced (avg)')
        ax1.plot(range(len(baseline_avg)), baseline_avg, color='blue', linewidth=3, label='Baseline (avg)')
        ax1.axhline(y=100, color='green', linestyle='--', alpha=0.7, label='Initial Balance')
        ax1.axhline(y=0, color='black', linestyle='--', alpha=0.7, label='Bankruptcy')
        
        ax1.set_title('Balance Progression Comparison')
        ax1.set_xlabel('Episode')
        ax1.set_ylabel('Balance (Gold)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Final balance distribution
        enhanced_finals = [r['final_balance'] for r in enhanced_results]
        baseline_finals = [r['final_balance'] for r in baseline_results]
        
        ax2.hist(baseline_finals, bins=15, alpha=0.6, label='Baseline', color='blue', density=True)
        ax2.hist(enhanced_finals, bins=15, alpha=0.6, label='Enhanced', color='red', density=True)
        ax2.axvline(x=100, color='green', linestyle='--', alpha=0.7, label='Initial Balance')
        ax2.set_title('Final Balance Distribution')
        ax2.set_xlabel('Final Balance (Gold)')
        ax2.set_ylabel('Density')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Performance metrics comparison
        metrics = ['Mean Final', 'Profitable %', 'Bankrupt %']
        enhanced_values = [
            comparison_stats['enhanced_mean'],
            comparison_stats['enhanced_profitable'] / comparison_stats['num_runs'] * 100,
            comparison_stats['enhanced_bankrupt'] / comparison_stats['num_runs'] * 100
        ]
        baseline_values = [
            comparison_stats['baseline_mean'],
            comparison_stats['baseline_profitable'] / comparison_stats['num_runs'] * 100,
            comparison_stats['baseline_bankrupt'] / comparison_stats['num_runs'] * 100
        ]
        
        x = np.arange(len(metrics))
        width = 0.35
        
        ax3.bar(x - width/2, baseline_values, width, label='Baseline', color='blue', alpha=0.7)
        ax3.bar(x + width/2, enhanced_values, width, label='Enhanced', color='red', alpha=0.7)
        ax3.set_title('Performance Metrics Comparison')
        ax3.set_xlabel('Metrics')
        ax3.set_ylabel('Value')
        ax3.set_xticks(x)
        ax3.set_xticklabels(metrics)
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')
        
        # Add value labels on bars
        for i, (b_val, e_val) in enumerate(zip(baseline_values, enhanced_values)):
            ax3.text(i - width/2, b_val + 1, f'{b_val:.1f}', ha='center', va='bottom')
            ax3.text(i + width/2, e_val + 1, f'{e_val:.1f}', ha='center', va='bottom')
        
        # Improvement summary
        improvement_text = f"""Comparison Summary:
Enhanced vs Baseline

Final Balance:
Enhanced: {comparison_stats['enhanced_mean']:.1f} ± {comparison_stats['enhanced_std']:.1f}
Baseline: {comparison_stats['baseline_mean']:.1f} ± {comparison_stats['baseline_std']:.1f}

Improvement: {comparison_stats['improvement_percent']:+.1f}%

Win Rate:
Enhanced: {comparison_stats['enhanced_profitable']}/{comparison_stats['num_runs']} ({comparison_stats['enhanced_profitable']/comparison_stats['num_runs']*100:.1f}%)
Baseline: {comparison_stats['baseline_profitable']}/{comparison_stats['num_runs']} ({comparison_stats['baseline_profitable']/comparison_stats['num_runs']*100:.1f}%)"""
        
        ax4.text(0.05, 0.95, improvement_text, transform=ax4.transAxes, fontsize=10,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax4.set_xlim(0, 1)
        ax4.set_ylim(0, 1)
        ax4.axis('off')
        ax4.set_title('Summary Statistics')
        
        plt.tight_layout()
        plt.savefig('enhanced_vs_baseline_simulation.png', dpi=300, bbox_inches='tight')
        plt.close()
        print("📊 Agent comparison plot saved to 'enhanced_vs_baseline_simulation.png'")
    
    def _calculate_average_curve(self, results: List[Dict], max_episodes: int) -> List[float]:
        """Calculate average balance curve across all runs."""
        padded_histories = []
        for result in results:
            history = result['balance_history']
            # Pad with final balance to make all curves same length
            padded = history + [history[-1]] * (max_episodes - len(history))
            padded_histories.append(padded)
        
        return np.mean(padded_histories, axis=0)

def main():
    """Main function to compare enhanced vs baseline agents."""
    print("🎯 Enhanced vs Baseline Agent Simulation Comparison")
    print("=" * 60)
    
    # Initialize simulator
    simulator = EnhancedGameplaySimulator()
    
    # Load both models
    if not simulator.load_models():
        print("\n❌ Could not load required models. Please ensure both models are trained:")
        print("   - Run 'python main.py' to train baseline model")
        print("   - Run 'python train_enhanced.py' to train enhanced model")
        return
    
    print("\n🎰 Running comprehensive agent comparison...")
    print("This will simulate realistic gameplay for both agents and compare performance.")
    print()
    
    # Run comparison
    comparison_stats = simulator.compare_agents(
        num_runs=25,
        initial_balance=100,
        max_episodes=500
    )
    
    # Save detailed results
    with open("agent_comparison_results.txt", 'w') as f:
        f.write("ENHANCED VS BASELINE AGENT COMPARISON\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("SIMULATION PARAMETERS:\n")
        f.write(f"Number of runs per agent: {comparison_stats['num_runs']}\n")
        f.write(f"Initial balance per run: 100 Gold\n")
        f.write(f"Maximum episodes per run: 500\n\n")
        
        f.write("ENHANCED AGENT RESULTS:\n")
        f.write(f"Mean final balance: {comparison_stats['enhanced_mean']:.2f} ± {comparison_stats['enhanced_std']:.2f} Gold\n")
        f.write(f"Profitable runs: {comparison_stats['enhanced_profitable']}/{comparison_stats['num_runs']} ({comparison_stats['enhanced_profitable']/comparison_stats['num_runs']*100:.1f}%)\n")
        f.write(f"Bankrupt runs: {comparison_stats['enhanced_bankrupt']}/{comparison_stats['num_runs']} ({comparison_stats['enhanced_bankrupt']/comparison_stats['num_runs']*100:.1f}%)\n\n")
        
        f.write("BASELINE AGENT RESULTS:\n")
        f.write(f"Mean final balance: {comparison_stats['baseline_mean']:.2f} ± {comparison_stats['baseline_std']:.2f} Gold\n")
        f.write(f"Profitable runs: {comparison_stats['baseline_profitable']}/{comparison_stats['num_runs']} ({comparison_stats['baseline_profitable']/comparison_stats['num_runs']*100:.1f}%)\n")
        f.write(f"Bankrupt runs: {comparison_stats['baseline_bankrupt']}/{comparison_stats['num_runs']} ({comparison_stats['baseline_bankrupt']/comparison_stats['num_runs']*100:.1f}%)\n\n")
        
        f.write("PERFORMANCE IMPROVEMENT:\n")
        f.write(f"Improvement: {comparison_stats['improvement_percent']:+.2f}%\n")
    
    print(f"\n🎉 Comparison complete! Generated files:")
    print(f"   📊 enhanced_vs_baseline_simulation.png - Comprehensive comparison visualization")
    print(f"   📄 agent_comparison_results.txt - Detailed comparison statistics")
    
    if comparison_stats['improvement_percent'] > 0:
        print(f"\n🚀 The enhanced agent shows {comparison_stats['improvement_percent']:.2f}% better performance!")
        print("The improvements from Double Q-Learning, adaptive exploration, and reward shaping are working!")
    else:
        print(f"\n🤔 The enhanced agent performed {abs(comparison_stats['improvement_percent']):.2f}% worse.")
        print("Consider further hyperparameter tuning or additional enhancement techniques.")

if __name__ == "__main__":
    main()