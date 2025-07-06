#!/usr/bin/env python3

"""
Train the enhanced Q-learning agent and compare with baseline.
"""

import numpy as np
import matplotlib.pyplot as plt
from game_rules import SlotMachine
from enhanced_q_learning import EnhancedQLearningAgent
from q_learning import QLearningAgent
from utils import plot_learning_curve, analyze_policy, save_results

def train_enhanced_agent():
    """Train the enhanced Q-learning agent."""
    print("🚀 Training Enhanced Q-Learning Agent")
    print("=" * 50)
    
    # Initialize the slot machine game
    game = SlotMachine()
    print(f"Game initialized:")
    print(f"  Symbols: {game.symbols}")
    print(f"  Abstract state space size: {game.get_abstract_state_space_size()}")
    print(f"  Action space size: {game.get_action_space_size()}")
    print()
    
    # Initialize enhanced Q-learning agent
    enhanced_agent = EnhancedQLearningAgent(
        state_space_size=game.get_abstract_state_space_size(),
        action_space_size=game.get_action_space_size(),
        learning_rate=0.25,
        discount_factor=0.99,
        epsilon=1.0,
        epsilon_decay=0.9995,
        epsilon_min=0.01,
        lr_decay=0.9999,
        lr_min=0.01,
        use_double_q=True,      # Enable Double Q-Learning
        adaptive_epsilon=True   # Enable adaptive exploration
    )
    
    print("Enhanced Q-learning agent initialized:")
    print(f"  Double Q-Learning: {enhanced_agent.use_double_q}")
    print(f"  Adaptive Epsilon: {enhanced_agent.adaptive_epsilon}")
    print(f"  Learning rate: {enhanced_agent.learning_rate}")
    print(f"  Discount factor: {enhanced_agent.discount_factor}")
    print(f"  Initial epsilon: {enhanced_agent.epsilon}")
    print()
    
    # Training
    num_episodes = 100000
    print(f"Starting enhanced training for {num_episodes} episodes...")
    
    enhanced_rewards = enhanced_agent.train(
        game=game,
        num_episodes=num_episodes,
        max_steps_per_episode=20,
        verbose=True,
        report_interval=5000
    )
    
    print("\nEnhanced training completed!")
    print("-" * 50)
    
    # Analyze results
    print("\nEnhanced Agent Analysis:")
    enhanced_policy = enhanced_agent.get_policy()
    analyze_policy(enhanced_policy, game, game.actions)
    
    # Evaluate final policy
    print("\nEnhanced Policy Evaluation:")
    mean_reward, std_reward = enhanced_agent.evaluate_policy(game, num_episodes=10000)
    print(f"Average reward per episode: {mean_reward:.4f} ± {std_reward:.4f}")
    
    # Save enhanced model and results
    enhanced_agent.save_model("enhanced_slot_machine_model.npz")
    save_results(enhanced_agent.q_table_a if enhanced_agent.use_double_q else enhanced_agent.q_table, 
                enhanced_rewards, enhanced_policy, "enhanced_slot_machine_results")
    
    return enhanced_rewards, enhanced_policy, mean_reward, std_reward

def compare_with_baseline():
    """Compare enhanced agent with baseline agent."""
    print("\n" + "="*60)
    print("📊 COMPARING ENHANCED VS BASELINE AGENT")
    print("="*60)
    
    game = SlotMachine()
    
    # Train baseline agent for comparison
    print("\n🔄 Training baseline agent for comparison...")
    baseline_agent = QLearningAgent(
        state_space_size=game.get_abstract_state_space_size(),
        action_space_size=game.get_action_space_size(),
        learning_rate=0.25,
        discount_factor=0.99,
        epsilon=1.0,
        epsilon_decay=0.9995,
        epsilon_min=0.01
    )
    
    baseline_rewards = baseline_agent.train(
        game=game,
        num_episodes=100000,
        verbose=False,  # Suppress output for comparison
        report_interval=10000
    )
    
    # Evaluate baseline
    baseline_mean, baseline_std = baseline_agent.evaluate_policy(game, num_episodes=10000)
    
    # Load enhanced results
    try:
        data = np.load("enhanced_slot_machine_results.npz")
        enhanced_rewards = data['rewards'].tolist()
    except FileNotFoundError:
        print("Enhanced results not found. Please run enhanced training first.")
        return
    
    # Train enhanced agent
    enhanced_rewards, enhanced_policy, enhanced_mean, enhanced_std = train_enhanced_agent()
    
    # Performance comparison
    print(f"\n📈 PERFORMANCE COMPARISON:")
    print(f"Baseline Agent:")
    print(f"  Final avg reward: {baseline_mean:.4f} ± {baseline_std:.4f}")
    print(f"  Best episode reward: {max(baseline_rewards):.2f}")
    print(f"  Final 1000 episodes avg: {np.mean(baseline_rewards[-1000:]):.4f}")
    
    print(f"\nEnhanced Agent:")
    print(f"  Final avg reward: {enhanced_mean:.4f} ± {enhanced_std:.4f}")
    print(f"  Best episode reward: {max(enhanced_rewards):.2f}")
    print(f"  Final 1000 episodes avg: {np.mean(enhanced_rewards[-1000:]):.4f}")
    
    # Calculate improvement
    improvement = (enhanced_mean - baseline_mean) / abs(baseline_mean) * 100
    print(f"\n🎯 IMPROVEMENT: {improvement:+.2f}%")
    
    # Create comparison plot
    create_comparison_plot(baseline_rewards, enhanced_rewards)
    
    return {
        'baseline_mean': baseline_mean,
        'baseline_std': baseline_std,
        'enhanced_mean': enhanced_mean,
        'enhanced_std': enhanced_std,
        'improvement_percent': improvement
    }

def create_comparison_plot(baseline_rewards, enhanced_rewards):
    """Create visualization comparing baseline vs enhanced learning curves."""
    plt.figure(figsize=(15, 10))
    
    # Learning curves comparison
    plt.subplot(2, 2, 1)
    window_size = 1000
    
    # Baseline moving average
    if len(baseline_rewards) >= window_size:
        baseline_ma = np.convolve(baseline_rewards, np.ones(window_size)/window_size, mode='valid')
        plt.plot(baseline_ma, label='Baseline Agent', color='blue', alpha=0.8)
    
    # Enhanced moving average
    if len(enhanced_rewards) >= window_size:
        enhanced_ma = np.convolve(enhanced_rewards, np.ones(window_size)/window_size, mode='valid')
        plt.plot(enhanced_ma, label='Enhanced Agent', color='red', alpha=0.8)
    
    plt.title('Learning Curves Comparison (Moving Average)')
    plt.xlabel('Episode')
    plt.ylabel('Average Reward')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Final performance comparison
    plt.subplot(2, 2, 2)
    final_baseline = np.mean(baseline_rewards[-1000:])
    final_enhanced = np.mean(enhanced_rewards[-1000:])
    
    agents = ['Baseline', 'Enhanced']
    performances = [final_baseline, final_enhanced]
    colors = ['blue', 'red']
    
    bars = plt.bar(agents, performances, color=colors, alpha=0.7)
    plt.title('Final Performance Comparison')
    plt.ylabel('Average Reward (Last 1000 Episodes)')
    
    # Add value labels on bars
    for bar, perf in zip(bars, performances):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f'{perf:.3f}', ha='center', va='bottom')
    
    plt.grid(True, alpha=0.3, axis='y')
    
    # Reward distribution comparison
    plt.subplot(2, 2, 3)
    plt.hist(baseline_rewards[-5000:], bins=50, alpha=0.6, label='Baseline', color='blue', density=True)
    plt.hist(enhanced_rewards[-5000:], bins=50, alpha=0.6, label='Enhanced', color='red', density=True)
    plt.title('Reward Distribution (Last 5000 Episodes)')
    plt.xlabel('Episode Reward')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Convergence comparison
    plt.subplot(2, 2, 4)
    # Calculate convergence: running average of last 100 episodes
    baseline_convergence = [np.mean(baseline_rewards[max(0, i-100):i+1]) for i in range(99, len(baseline_rewards))]
    enhanced_convergence = [np.mean(enhanced_rewards[max(0, i-100):i+1]) for i in range(99, len(enhanced_rewards))]
    
    plt.plot(baseline_convergence, label='Baseline', color='blue', alpha=0.8)
    plt.plot(enhanced_convergence, label='Enhanced', color='red', alpha=0.8)
    plt.title('Convergence Comparison (Running Average)')
    plt.xlabel('Episode')
    plt.ylabel('Running Average Reward')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('enhanced_vs_baseline_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("📊 Comparison plot saved to 'enhanced_vs_baseline_comparison.png'")

def main():
    """Main function to train enhanced agent and compare performance."""
    print("🎯 Enhanced Q-Learning Agent Training & Comparison")
    print("=" * 60)
    
    # Train enhanced agent
    enhanced_rewards, enhanced_policy, enhanced_mean, enhanced_std = train_enhanced_agent()
    
    # Generate learning curve visualization
    plot_learning_curve(enhanced_rewards, window_size=1000, filename="enhanced_learning_curve.png")
    
    print("\n✅ Enhanced agent training complete!")
    print("🔍 To compare with baseline agent, the baseline will be trained automatically.")
    print("📊 Generating comparison analysis...")
    
    # Compare with baseline
    comparison_results = compare_with_baseline()
    
    print(f"\n🎉 Analysis complete! Generated files:")
    print(f"   📊 enhanced_learning_curve.png - Enhanced agent learning curve")
    print(f"   📊 enhanced_vs_baseline_comparison.png - Performance comparison")
    print(f"   💾 enhanced_slot_machine_model.npz - Trained enhanced model")
    print(f"   📄 enhanced_slot_machine_results.npz - Training data")
    
    if comparison_results['improvement_percent'] > 0:
        print(f"\n🚀 The enhanced agent achieved {comparison_results['improvement_percent']:.2f}% improvement!")
    else:
        print(f"\n⚠️  The enhanced agent performed {abs(comparison_results['improvement_percent']):.2f}% worse.")
        print("Consider adjusting hyperparameters or trying different enhancement techniques.")

if __name__ == "__main__":
    main()