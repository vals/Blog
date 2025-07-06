#!/usr/bin/env python3
"""
Test script for the Hybrid Expected Value Strategy.
Compares hybrid EV agent with naive cash-out strategy and standard Q-learning.
"""

import numpy as np
import json
from game_rules import SlotMachine
from hybrid_ev_agent import HybridEVAgent
from q_learning import QLearningAgent

def test_hybrid_strategy():
    """Test the hybrid expected value strategy."""
    print("="*60)
    print("HYBRID EXPECTED VALUE STRATEGY TEST")
    print("="*60)
    
    # Initialize game
    game = SlotMachine()
    
    # Test 1: Basic Hybrid Strategy Performance
    print("\n1. Testing Basic Hybrid Strategy Performance")
    print("-" * 40)
    
    hybrid_agent = HybridEVAgent(game, ev_threshold=0.0)
    hybrid_results = hybrid_agent.evaluate_strategy(num_episodes=10000)
    
    print(f"Hybrid Strategy Results:")
    print(f"  Average Reward: {hybrid_results['avg_reward']:.3f}")
    print(f"  Standard Deviation: {hybrid_results['std_reward']:.3f}")
    print(f"  Profitable Episodes: {hybrid_results['profitable_episodes']}/{10000}")
    print(f"  Profitable Rate: {hybrid_results['profitable_rate']:.1f}%")
    print(f"  Best Reward: {hybrid_results['best_reward']:.0f}")
    print(f"  Worst Reward: {hybrid_results['worst_reward']:.0f}")
    print(f"  Cash Out Percentage: {hybrid_results['cash_out_percentage']:.1f}%")
    print(f"  EV Override Percentage: {hybrid_results['ev_override_percentage']:.1f}%")
    
    # Test 2: Compare with Naive Strategy
    print("\n2. Comparing with Naive Cash-Out Strategy")
    print("-" * 40)
    
    comparison = hybrid_agent.compare_with_naive_strategy(num_episodes=10000)
    
    print(f"Strategy Comparison:")
    print(f"  Hybrid Avg Reward: {comparison['hybrid_strategy']['avg_reward']:.3f}")
    print(f"  Naive Avg Reward: {comparison['naive_strategy']['avg_reward']:.3f}")
    print(f"  Improvement: {comparison['performance_difference']['reward_improvement']:.3f}")
    print(f"  Improvement %: {comparison['performance_difference']['improvement_percentage']:.1f}%")
    print(f"  Profitable Rate Improvement: {comparison['performance_difference']['profitable_rate_improvement']:.1f}%")
    
    # Test 3: Analyze EV Override Decisions
    print("\n3. Analyzing EV Override Decisions")
    print("-" * 40)
    
    override_analysis = hybrid_agent.analyze_ev_overrides()
    
    if 'message' not in override_analysis:
        print(f"EV Override Analysis:")
        print(f"  Total Overrides: {override_analysis['total_overrides']}")
        print(f"  Average EV Advantage: {override_analysis['avg_ev_advantage']:.3f}")
        print(f"  Max EV Advantage: {override_analysis['max_ev_advantage']:.3f}")
        print(f"  Average Cash-Out Value: {override_analysis['avg_cash_out_value']:.3f}")
        print(f"  Average Chosen EV: {override_analysis['avg_chosen_ev']:.3f}")
        
        print(f"\n  Most Common Override Actions:")
        for action, count in list(override_analysis['action_distribution'].items())[:5]:
            print(f"    {action}: {count} times")
        
        print(f"\n  Most Common State Patterns:")
        for pattern, count in list(override_analysis['common_state_patterns'].items())[:5]:
            print(f"    {pattern}: {count} times")
    else:
        print(f"  {override_analysis['message']}")
    
    # Test 4: Different EV Thresholds
    print("\n4. Testing Different EV Thresholds")
    print("-" * 40)
    
    thresholds = [0.0, 0.5, 1.0, 2.0]
    threshold_results = {}
    
    for threshold in thresholds:
        print(f"\nTesting EV Threshold: {threshold}")
        agent = HybridEVAgent(game, ev_threshold=threshold)
        results = agent.evaluate_strategy(num_episodes=5000)
        threshold_results[threshold] = results
        
        print(f"  Avg Reward: {results['avg_reward']:.3f}")
        print(f"  Profitable Rate: {results['profitable_rate']:.1f}%")
        print(f"  EV Override Rate: {results['ev_override_percentage']:.1f}%")
    
    # Test 5: Compare with Q-Learning (if model exists)
    print("\n5. Comparing with Q-Learning Agent")
    print("-" * 40)
    
    try:
        # Try to load existing Q-learning model
        q_agent = QLearningAgent(
            state_space_size=game.get_abstract_state_space_size(),
            action_space_size=game.get_action_space_size()
        )
        
        try:
            q_agent.load_model('slot_machine_model.npz')
            print("Loaded existing Q-learning model")
            
            # Evaluate Q-learning performance
            q_avg_reward, q_std_reward = q_agent.evaluate_policy(game, num_episodes=5000)
            
            print(f"Q-Learning Results:")
            print(f"  Average Reward: {q_avg_reward:.3f}")
            print(f"  Standard Deviation: {q_std_reward:.3f}")
            
            print(f"\nPerformance Comparison:")
            print(f"  Hybrid Strategy: {hybrid_results['avg_reward']:.3f}")
            print(f"  Q-Learning: {q_avg_reward:.3f}")
            print(f"  Difference: {hybrid_results['avg_reward'] - q_avg_reward:.3f}")
            
        except FileNotFoundError:
            print("No existing Q-learning model found. Training new model...")
            
            # Train a new Q-learning model for comparison
            q_rewards = q_agent.train(game, num_episodes=50000, verbose=False)
            q_avg_reward, q_std_reward = q_agent.evaluate_policy(game, num_episodes=5000)
            
            print(f"Q-Learning Results (after training):")
            print(f"  Average Reward: {q_avg_reward:.3f}")
            print(f"  Standard Deviation: {q_std_reward:.3f}")
            
            print(f"\nPerformance Comparison:")
            print(f"  Hybrid Strategy: {hybrid_results['avg_reward']:.3f}")
            print(f"  Q-Learning: {q_avg_reward:.3f}")
            print(f"  Difference: {hybrid_results['avg_reward'] - q_avg_reward:.3f}")
            
    except Exception as e:
        print(f"Could not compare with Q-learning: {e}")
    
    # Save results
    print("\n6. Saving Results")
    print("-" * 40)
    
    # Save hybrid agent results
    hybrid_agent.save_results('hybrid_ev_strategy_results.npz')
    
    # Save comprehensive comparison
    final_results = {
        'hybrid_strategy': hybrid_results,
        'naive_comparison': comparison,
        'ev_override_analysis': override_analysis,
        'threshold_comparison': threshold_results
    }
    
    with open('hybrid_strategy_analysis.json', 'w') as f:
        # Convert numpy types to Python types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            return obj
        
        json_results = json.loads(json.dumps(final_results, default=convert_numpy))
        json.dump(json_results, f, indent=2)
    
    print("Results saved to:")
    print("  - hybrid_ev_strategy_results.npz")
    print("  - hybrid_strategy_analysis.json")
    
    print("\n" + "="*60)
    print("HYBRID STRATEGY TEST COMPLETED")
    print("="*60)
    
    return final_results

if __name__ == "__main__":
    results = test_hybrid_strategy()
    
    # Quick summary
    print("\nQUICK SUMMARY:")
    print(f"Hybrid Strategy Average Reward: {results['hybrid_strategy']['avg_reward']:.3f}")
    print(f"Naive Strategy Average Reward: {results['naive_comparison']['naive_strategy']['avg_reward']:.3f}")
    print(f"Performance Improvement: {results['naive_comparison']['performance_difference']['reward_improvement']:.3f}")
    
    if results['hybrid_strategy']['ev_override_percentage'] > 0:
        print(f"EV Override Rate: {results['hybrid_strategy']['ev_override_percentage']:.1f}%")
        print("The hybrid strategy found positive EV opportunities!")
    else:
        print("No positive EV opportunities found - strategy defaulted to cash-out only.")