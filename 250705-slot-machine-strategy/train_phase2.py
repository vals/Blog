#!/usr/bin/env python3

"""
Train Phase 2 advanced Q-learning agents and compare with Phase 1.
Includes UCB exploration, Prioritized Experience Replay, Multi-step learning, and Ensemble methods.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Dict
from game_rules import SlotMachine
from advanced_q_learning import AdvancedQLearningAgent
from ensemble_q_learning import EnsembleQLearningAgent
from enhanced_q_learning import EnhancedQLearningAgent
from q_learning import QLearningAgent
from utils import analyze_policy, save_results

def train_advanced_agent():
    """Train the advanced Q-learning agent with all Phase 2 features."""
    print("🚀 Training Advanced Q-Learning Agent (Phase 2)")
    print("=" * 60)
    
    game = SlotMachine()
    print(f"Game initialized:")
    print(f"  Abstract state space size: {game.get_abstract_state_space_size()}")
    print(f"  Action space size: {game.get_action_space_size()}")
    print()
    
    # Initialize advanced agent
    advanced_agent = AdvancedQLearningAgent(
        state_space_size=game.get_abstract_state_space_size(),
        action_space_size=game.get_action_space_size(),
        learning_rate=0.25,
        discount_factor=0.99,
        epsilon=1.0,
        epsilon_decay=0.9995,
        epsilon_min=0.01,
        lr_decay=0.9999,
        lr_min=0.01,
        use_double_q=True,
        use_ucb=True,
        ucb_c=2.0,
        use_replay=True,
        replay_capacity=50000,
        batch_size=32,
        n_step=3
    )
    
    print("Advanced Q-learning agent initialized:")
    print(f"  Double Q-Learning: {advanced_agent.use_double_q}")
    print(f"  UCB Exploration: {advanced_agent.use_ucb} (c={advanced_agent.ucb_c})")
    if advanced_agent.use_replay:
        print(f"  Prioritized Replay: {advanced_agent.use_replay} (capacity={advanced_agent.replay_buffer.capacity})")
    else:
        print(f"  Prioritized Replay: {advanced_agent.use_replay}")
    print(f"  N-step Learning: {advanced_agent.n_step}")
    print(f"  Learning rate: {advanced_agent.learning_rate}")
    print(f"  Discount factor: {advanced_agent.discount_factor}")
    print()
    
    # Training
    num_episodes = 100000
    print(f"Starting advanced training for {num_episodes} episodes...")
    
    advanced_rewards = advanced_agent.train(
        game=game,
        num_episodes=num_episodes,
        max_steps_per_episode=20,
        verbose=True,
        report_interval=5000
    )
    
    print("\nAdvanced training completed!")
    print("-" * 50)
    
    # Analyze results
    print("\nAdvanced Agent Analysis:")
    advanced_policy = advanced_agent.get_policy()
    analyze_policy(advanced_policy, game, game.actions)
    
    # Evaluate final policy
    print("\nAdvanced Policy Evaluation:")
    mean_reward, std_reward = advanced_agent.evaluate_policy(game, num_episodes=10000)
    print(f"Average reward per episode: {mean_reward:.4f} ± {std_reward:.4f}")
    
    # Save model and results
    advanced_agent.save_model("advanced_slot_machine_model.npz")
    save_results(advanced_agent.q_table_a if advanced_agent.use_double_q else advanced_agent.q_table, 
                advanced_rewards, advanced_policy, "advanced_slot_machine_results")
    
    return advanced_rewards, advanced_policy, mean_reward, std_reward

def train_ensemble_agent():
    """Train the ensemble Q-learning agent."""
    print("\n🤖 Training Ensemble Q-Learning Agent")
    print("=" * 60)
    
    game = SlotMachine()
    
    # Initialize ensemble agent
    ensemble_agent = EnsembleQLearningAgent(
        state_space_size=game.get_abstract_state_space_size(),
        action_space_size=game.get_action_space_size()
    )
    
    print("Ensemble agent initialized with 5 diverse agents")
    print()
    
    # Training
    num_episodes = 100000
    print(f"Starting ensemble training for {num_episodes} episodes...")
    
    ensemble_rewards = ensemble_agent.train(
        game=game,
        num_episodes=num_episodes,
        max_steps_per_episode=20,
        verbose=True,
        report_interval=5000
    )
    
    print("\nEnsemble training completed!")
    print("-" * 50)
    
    # Analyze results
    print("\nEnsemble Agent Analysis:")
    ensemble_policy = ensemble_agent.get_policy()
    analyze_policy(ensemble_policy, game, game.actions)
    
    # Evaluate final policy
    print("\nEnsemble Policy Evaluation:")
    mean_reward, std_reward = ensemble_agent.evaluate_policy(game, num_episodes=10000)
    print(f"Average reward per episode: {mean_reward:.4f} ± {std_reward:.4f}")
    
    # Save model and results
    ensemble_agent.save_model("ensemble_slot_machine_model")
    save_results(np.zeros((game.get_abstract_state_space_size(), game.get_action_space_size())), 
                ensemble_rewards, ensemble_policy, "ensemble_slot_machine_results")
    
    return ensemble_rewards, ensemble_policy, mean_reward, std_reward

def compare_all_phases():
    """Compare all phases: Baseline, Phase 1 (Enhanced), Phase 2 (Advanced), and Ensemble."""
    print("\n📊 COMPREHENSIVE PHASE COMPARISON")
    print("=" * 60)
    
    game = SlotMachine()
    results = {}
    
    # Train all agents for fair comparison
    agents_to_train = [
        ("Baseline", QLearningAgent, {
            'learning_rate': 0.25,
            'discount_factor': 0.99,
            'epsilon': 1.0,
            'epsilon_decay': 0.9995,
            'epsilon_min': 0.01
        }),
        ("Enhanced (Phase 1)", EnhancedQLearningAgent, {
            'learning_rate': 0.25,
            'discount_factor': 0.99,
            'epsilon': 1.0,
            'epsilon_decay': 0.9995,
            'epsilon_min': 0.01,
            'use_double_q': True,
            'adaptive_epsilon': True
        })
    ]
    
    print("\n🔄 Training agents for comparison...")
    
    for name, agent_class, params in agents_to_train:
        print(f"\nTraining {name}...")
        
        agent = agent_class(
            state_space_size=game.get_abstract_state_space_size(),
            action_space_size=game.get_action_space_size(),
            **params
        )
        
        rewards = agent.train(
            game=game,
            num_episodes=100000,
            verbose=False,
            report_interval=20000
        )
        
        mean_reward, std_reward = agent.evaluate_policy(game, num_episodes=5000)
        
        results[name] = {
            'rewards': rewards,
            'mean_reward': mean_reward,
            'std_reward': std_reward,
            'final_1000': np.mean(rewards[-1000:]),
            'best_reward': max(rewards)
        }
        
        print(f"  {name} - Final performance: {mean_reward:.4f} ± {std_reward:.4f}")
    
    # Train Phase 2 agents
    print(f"\nTraining Advanced (Phase 2)...")
    adv_rewards, adv_policy, adv_mean, adv_std = train_advanced_agent()
    results["Advanced (Phase 2)"] = {
        'rewards': adv_rewards,
        'mean_reward': adv_mean,
        'std_reward': adv_std,
        'final_1000': np.mean(adv_rewards[-1000:]),
        'best_reward': max(adv_rewards)
    }
    
    print(f"\nTraining Ensemble...")
    ens_rewards, ens_policy, ens_mean, ens_std = train_ensemble_agent()
    results["Ensemble"] = {
        'rewards': ens_rewards,
        'mean_reward': ens_mean,
        'std_reward': ens_std,
        'final_1000': np.mean(ens_rewards[-1000:]),
        'best_reward': max(ens_rewards)
    }
    
    # Create comprehensive comparison
    create_phase_comparison_plot(results)
    
    # Print comparison table
    print("\n📈 FINAL PERFORMANCE COMPARISON:")
    print("-" * 80)
    print(f"{'Agent':<25} {'Final Eval':<15} {'Last 1000':<12} {'Best':<8} {'Improvement':<12}")
    print("-" * 80)
    
    baseline_mean = results["Baseline"]['mean_reward']
    
    for name, data in results.items():
        improvement = (data['mean_reward'] - baseline_mean) / abs(baseline_mean) * 100
        print(f"{name:<25} {data['mean_reward']:>6.3f} ± {data['std_reward']:<4.3f} "
              f"{data['final_1000']:>8.3f}    {data['best_reward']:>6.1f}  {improvement:>+6.1f}%")
    
    return results

def create_phase_comparison_plot(results: Dict):
    """Create comprehensive visualization comparing all phases."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    
    # Learning curves comparison
    window_size = 1000
    for i, (name, data) in enumerate(results.items()):
        rewards = data['rewards']
        if len(rewards) >= window_size:
            moving_avg = np.convolve(rewards, np.ones(window_size)/window_size, mode='valid')
            ax1.plot(moving_avg, label=name, color=colors[i], alpha=0.8)
    
    ax1.set_title('Learning Curves Comparison (Moving Average)', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Episode')
    ax1.set_ylabel('Average Reward')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Final performance comparison
    names = list(results.keys())
    performances = [results[name]['mean_reward'] for name in names]
    errors = [results[name]['std_reward'] for name in names]
    
    bars = ax2.bar(names, performances, yerr=errors, capsize=5, color=colors[:len(names)], alpha=0.7)
    ax2.set_title('Final Performance Comparison', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Average Reward (Evaluation)')
    ax2.tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar, perf in zip(bars, performances):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{perf:.3f}', ha='center', va='bottom')
    
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Improvement comparison
    baseline_perf = results["Baseline"]['mean_reward']
    improvements = [(results[name]['mean_reward'] - baseline_perf) / abs(baseline_perf) * 100 
                   for name in names]
    
    colors_improvement = ['gray' if imp < 0 else 'green' for imp in improvements]
    bars = ax3.bar(names, improvements, color=colors_improvement, alpha=0.7)
    ax3.set_title('Improvement Over Baseline (%)', fontsize=14, fontweight='bold')
    ax3.set_ylabel('Improvement (%)')
    ax3.tick_params(axis='x', rotation=45)
    ax3.axhline(y=0, color='black', linestyle='-', alpha=0.5)
    
    # Add value labels
    for bar, imp in zip(bars, improvements):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + (1 if imp >= 0 else -3),
                f'{imp:+.1f}%', ha='center', va='bottom' if imp >= 0 else 'top')
    
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Best episode rewards comparison
    best_rewards = [results[name]['best_reward'] for name in names]
    bars = ax4.bar(names, best_rewards, color=colors[:len(names)], alpha=0.7)
    ax4.set_title('Best Episode Reward Achieved', fontsize=14, fontweight='bold')
    ax4.set_ylabel('Best Reward')
    ax4.tick_params(axis='x', rotation=45)
    
    # Add value labels
    for bar, best in zip(bars, best_rewards):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{best:.1f}', ha='center', va='bottom')
    
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig('phase_comparison_comprehensive.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("📊 Comprehensive phase comparison saved to 'phase_comparison_comprehensive.png'")

def main():
    """Main function to train Phase 2 agents and compare all phases."""
    print("🎯 Phase 2 Advanced Q-Learning Training & Comprehensive Comparison")
    print("=" * 70)
    
    # Run comprehensive comparison
    results = compare_all_phases()
    
    # Save detailed comparison results
    with open("phase_comparison_results.txt", 'w') as f:
        f.write("COMPREHENSIVE PHASE COMPARISON RESULTS\n")
        f.write("=" * 50 + "\n\n")
        
        f.write("TRAINING PARAMETERS:\n")
        f.write("Episodes per agent: 100,000\n")
        f.write("Evaluation episodes: 5,000-10,000\n\n")
        
        baseline_mean = results["Baseline"]['mean_reward']
        
        for name, data in results.items():
            f.write(f"{name.upper()}:\n")
            f.write(f"  Final evaluation: {data['mean_reward']:.4f} ± {data['std_reward']:.4f}\n")
            f.write(f"  Last 1000 episodes: {data['final_1000']:.4f}\n")
            f.write(f"  Best episode: {data['best_reward']:.2f}\n")
            
            improvement = (data['mean_reward'] - baseline_mean) / abs(baseline_mean) * 100
            f.write(f"  Improvement over baseline: {improvement:+.2f}%\n\n")
    
    print(f"\n🎉 Phase 2 analysis complete! Generated files:")
    print(f"   📊 phase_comparison_comprehensive.png - Complete performance comparison")
    print(f"   📄 phase_comparison_results.txt - Detailed comparison statistics")
    print(f"   💾 advanced_slot_machine_model.npz - Advanced agent model")
    print(f"   💾 ensemble_slot_machine_model_*.npz - Ensemble agent models")
    
    # Summary
    best_agent = max(results.keys(), key=lambda k: results[k]['mean_reward'])
    best_improvement = (results[best_agent]['mean_reward'] - baseline_mean) / abs(baseline_mean) * 100
    
    print(f"\n🏆 BEST PERFORMING AGENT: {best_agent}")
    print(f"   Performance: {results[best_agent]['mean_reward']:.4f} ± {results[best_agent]['std_reward']:.4f}")
    print(f"   Improvement: {best_improvement:+.2f}% over baseline")
    
    if best_improvement > 10:
        print("\n🚀 Significant improvement achieved with advanced techniques!")
    elif best_improvement > 0:
        print("\n✅ Modest but positive improvement with advanced methods.")
    else:
        print("\n🤔 Advanced methods didn't outperform simpler approaches.")
        print("This might indicate the problem is well-suited to basic Q-learning,")
        print("or that hyperparameter tuning could further improve results.")

if __name__ == "__main__":
    main()