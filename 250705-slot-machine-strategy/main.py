#!/usr/bin/env python3

import numpy as np
from game_rules import SlotMachine
from q_learning import QLearningAgent
from utils import (
    print_q_table_summary, 
    get_policy_from_q_table, 
    plot_learning_curve,
    plot_policy_distribution,
    analyze_policy,
    save_results,
    load_results
)

def main():
    print("Slot Machine Q-Learning Strategy Optimization")
    print("=" * 50)
    
    # Initialize the slot machine game
    game = SlotMachine()
    print(f"Game initialized:")
    print(f"  Symbols: {game.symbols}")
    print(f"  Number of reels: {game.num_reels}")
    print(f"  Original state space size: {game.get_state_space_size()}")
    print(f"  Abstract state space size: {game.get_abstract_state_space_size()}")
    print(f"  Action space size: {game.get_action_space_size()}")
    print(f"  Actions: {game.actions}")
    print()
    
    # Initialize Q-learning agent with conservative parameters to prevent degradation
    agent = QLearningAgent(
        state_space_size=game.get_abstract_state_space_size(),  # Use abstract state space
        action_space_size=game.get_action_space_size(),
        learning_rate=0.15,        # More moderate learning rate
        discount_factor=0.95,      # Less aggressive future reward weighting
        epsilon=1.0,
        epsilon_decay=0.99992,     # Very conservative decay to maintain exploration
        epsilon_min=0.05           # Higher minimum to always maintain some exploration
    )
    
    print("Q-learning agent initialized with enhanced exploration:")
    print(f"  Learning rate: {agent.learning_rate} (with decay {agent.lr_decay})")
    print(f"  Discount factor: {agent.discount_factor}")
    print(f"  Initial epsilon: {agent.epsilon}")
    print(f"  Epsilon decay: {agent.epsilon_decay}")
    print(f"  Minimum epsilon: {agent.epsilon_min}")
    print(f"  Enhanced parameters for better convergence")
    print()
    
    # Training
    num_episodes = 100_000
    print(f"Starting training for {num_episodes} episodes...")
    
    rewards = agent.train(
        game=game,
        num_episodes=num_episodes,
        max_steps_per_episode=20,  # Up to 1 initial spin + 5 respins + decisions
        verbose=True,
        report_interval=5000
    )
    
    print("\nTraining completed!")
    print("-" * 50)
    
    # Analyze results
    print("\nQ-Table Analysis:")
    print_q_table_summary(agent.q_table, top_n=10)
    
    print("\nPolicy Analysis:")
    policy = agent.get_policy()
    analyze_policy(policy, game, game.actions)
    
    # Analyze strategy focus
    print("\nStrategy Focus Analysis:")
    focus_stats = agent.analyze_strategy_focus()
    print(f"  Mean Q-value std deviation: {focus_stats['mean_q_std']:.3f}")
    print(f"  Mean action entropy: {focus_stats['mean_action_entropy']:.3f}")
    print(f"  States with focused preferences: {focus_stats['focused_states_pct']:.1f}%")
    print(f"  Decisive states (low entropy): {focus_stats['decisive_states_pct']:.1f}%")
    
    # Evaluate final policy
    print("\nPolicy Evaluation:")
    mean_reward, std_reward = agent.evaluate_policy(game, num_episodes=10000)
    print(f"Average reward per spin: {mean_reward:.4f} ± {std_reward:.4f}")
    
    # Calculate average episode cost and return
    print(f"Average net reward per episode: {mean_reward:.4f} ± {std_reward:.4f}")
    print(f"Average episode includes initial spin (1 Gold) + potential respins (1 Gold each)")
    
    # Analyze action preferences
    action_distribution = {}
    for action in policy:
        action_distribution[action] = action_distribution.get(action, 0) + 1
    
    print(f"\nAction preferences in learned policy:")
    for action, count in action_distribution.items():
        percentage = (count / len(policy)) * 100
        action_desc = game.actions[action]['description']
        print(f"  {action_desc}: {percentage:.1f}% of states")
    
    # Generate visualizations
    print("\nGenerating visualizations...")
    try:
        plot_learning_curve(rewards, window_size=1000, filename="learning_curve.png")
        plot_policy_distribution(policy, game.actions, filename="policy_distribution.png")
        print("Visualizations saved as PNG files")
    except ImportError as e:
        print(f"Matplotlib not available - skipping visualizations: {e}")
    
    # Save results
    save_results(agent.q_table, rewards, policy, "slot_machine_results")
    agent.save_model("slot_machine_model.npz")
    
    print("\nTraining and analysis complete!")
    print("Results saved to 'slot_machine_results.npz' and 'slot_machine_model.npz'")
    print("Visualizations saved as 'learning_curve.png' and 'policy_distribution.png'")

def demo_game():
    print("\nDemo: Playing the slot machine")
    print("-" * 30)
    
    game = SlotMachine()
    
    for i in range(5):
        state = game.reset()
        print(f"Spin {i+1}: {state}")
        
        for action in range(game.get_action_space_size()):
            reward, next_state, _ = game.step(action)
            bet = game.actions[action]['bet']
            print(f"  Action {action} (bet {bet}): Reward = {reward}, Next state = {next_state}")
            game.current_state = state  # Reset to original state for fair comparison

def analyze_saved_results():
    try:
        q_table, rewards, policy = load_results("slot_machine_results.npz")
        print("\nLoaded saved results:")
        print_q_table_summary(q_table)
        plot_learning_curve(rewards)
    except FileNotFoundError:
        print("No saved results found. Run main() first.")

if __name__ == "__main__":
    # Run the main training
    main()
    
    # Uncomment to run demo
    # demo_game()
    
    # Uncomment to analyze saved results
    # analyze_saved_results()