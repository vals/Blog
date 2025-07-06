#!/usr/bin/env python3
"""
Test script for the Optimal Expected Value Q-Table.
Creates a theoretically optimal Q-table and tests its performance.
"""

import numpy as np
import json
from game_rules import SlotMachine
from optimal_ev_qtable import OptimalEVQTable, OptimalEVAgent

def test_optimal_qtable():
    """Test the optimal Q-table implementation."""
    print("="*60)
    print("OPTIMAL EXPECTED VALUE Q-TABLE TEST")
    print("="*60)
    
    # Initialize game
    game = SlotMachine()
    
    # Test 1: Build Optimal Q-Table
    print("\n1. Building Optimal Q-Table")
    print("-" * 40)
    
    qtable_builder = OptimalEVQTable(game)
    
    # Build the Q-table (this will take some time)
    q_table = qtable_builder.build_optimal_qtable()
    
    print(f"Q-table shape: {q_table.shape}")
    print(f"Q-table built successfully!")
    
    # Test 2: Analyze Q-Table
    print("\n2. Analyzing Q-Table")
    print("-" * 40)
    
    analysis = qtable_builder.analyze_qtable()
    print(f"Analysis Results:")
    print(f"  Total states: {analysis['total_states']}")
    print(f"  Total actions: {analysis['total_actions']}")
    print(f"  Positive EV actions: {analysis['positive_ev_actions']}")
    print(f"  Positive EV percentage: {analysis['positive_ev_percentage']:.1f}%")
    print(f"  States with positive EV respins: {analysis['states_with_positive_ev_respins']}")
    print(f"  Percentage of states with positive EV: {analysis['percentage_states_with_positive_ev']:.1f}%")
    print(f"  Value range: {analysis['value_statistics']['min_value']:.1f} to {analysis['value_statistics']['max_value']:.1f}")
    print(f"  Mean value: {analysis['value_statistics']['mean_value']:.3f}")
    
    # Test 3: Save Q-Table
    print("\n3. Saving Q-Table")
    print("-" * 40)
    
    qtable_builder.save_qtable('optimal_qtable.npz')
    
    # Test 4: Test Optimal Agent
    print("\n4. Testing Optimal Agent Performance")
    print("-" * 40)
    
    optimal_agent = OptimalEVAgent(game, qtable_builder)
    optimal_results = optimal_agent.evaluate_strategy(num_episodes=10000)
    
    print(f"Optimal Agent Results:")
    print(f"  Average Reward: {optimal_results['avg_reward']:.3f}")
    print(f"  Standard Deviation: {optimal_results['std_reward']:.3f}")
    print(f"  Profitable Episodes: {optimal_results['profitable_episodes']}/10000")
    print(f"  Profitable Rate: {optimal_results['profitable_rate']:.1f}%")
    print(f"  Best Reward: {optimal_results['best_reward']:.0f}")
    print(f"  Worst Reward: {optimal_results['worst_reward']:.0f}")
    print(f"  Cash Out Percentage: {optimal_results['cash_out_percentage']:.1f}%")
    print(f"  Respin Percentage: {optimal_results['respin_percentage']:.1f}%")
    
    # Test 5: Compare with Naive Strategy
    print("\n5. Comparing with Naive Strategy")
    print("-" * 40)
    
    # Run naive strategy for comparison
    naive_rewards = []
    print("Running Naive Strategy...")
    for _ in range(10000):
        game.reset()
        reward, _, _ = game.step(8)  # Always cash out
        naive_rewards.append(reward)
    
    naive_avg = np.mean(naive_rewards)
    naive_std = np.std(naive_rewards)
    naive_profitable = np.sum(np.array(naive_rewards) > 0)
    
    print(f"Strategy Comparison:")
    print(f"  Optimal Strategy: {optimal_results['avg_reward']:.3f}")
    print(f"  Naive Strategy: {naive_avg:.3f}")
    print(f"  Improvement: {optimal_results['avg_reward'] - naive_avg:.3f}")
    print(f"  Improvement %: {((optimal_results['avg_reward'] - naive_avg) / abs(naive_avg)) * 100:.1f}%")
    
    # Test 6: Find Best States
    print("\n6. Finding States with Highest Positive EV")
    print("-" * 40)
    
    # Find some high-value states
    high_value_states = []
    for state_idx in range(q_table.shape[0]):
        cash_out_val = q_table[state_idx, 8]
        respin_vals = q_table[state_idx, :8]
        valid_respins = respin_vals[respin_vals > float('-inf')]
        
        if len(valid_respins) > 0:
            best_respin = np.max(valid_respins)
            if best_respin > cash_out_val:
                ev_advantage = best_respin - cash_out_val
                high_value_states.append((state_idx, cash_out_val, best_respin, ev_advantage))
    
    # Sort by EV advantage
    high_value_states.sort(key=lambda x: x[3], reverse=True)
    
    print(f"Top 10 States with Highest EV Advantage:")
    for i, (state_idx, cash_out, best_respin, advantage) in enumerate(high_value_states[:10]):
        print(f"  {i+1}. State {state_idx}: Cash out={cash_out:.1f}, Best respin={best_respin:.1f}, Advantage={advantage:.1f}")
    
    # Test 7: Manual State Examples
    print("\n7. Testing Specific High-Value States")
    print("-" * 40)
    
    # Test some specific scenarios
    test_scenarios = [
        # [Blank, Coin, Coin Stack, Snake, Net, x2, Clover, Crown]
        ([0, 0, 0, 0, 0, 0, 0, 4], 0, "4 Crowns"),  # Jackpot
        ([0, 0, 0, 0, 0, 0, 0, 3], 0, "3 Crowns"),  # Near jackpot
        ([0, 0, 0, 0, 0, 0, 4, 0], 0, "4 Clovers"),  # All clovers
        ([0, 0, 0, 0, 0, 0, 3, 0], 0, "3 Clovers"),  # Near all clovers
        ([0, 4, 0, 0, 0, 0, 0, 0], 0, "4 Coins"),    # All coins
        ([0, 3, 0, 0, 0, 0, 0, 0], 0, "3 Coins"),    # Near all coins
        ([0, 0, 4, 0, 0, 0, 0, 0], 0, "4 Stacks"),   # All stacks
        ([0, 0, 3, 0, 0, 0, 0, 0], 0, "3 Stacks"),   # Near all stacks
        ([0, 0, 0, 2, 2, 0, 0, 0], 0, "2 Snakes + 2 Nets"),  # Snake-net combo
    ]
    
    for symbol_counts, respins, description in test_scenarios:
        try:
            action, value = qtable_builder.get_optimal_action(symbol_counts, respins)
            cash_out = qtable_builder._calculate_payout_from_counts(symbol_counts)
            action_name = "Cash Out" if action == 8 else game.actions[action]['description']
            
            print(f"  {description}:")
            print(f"    Cash out value: {cash_out:.1f}")
            print(f"    Optimal action: {action_name}")
            print(f"    Expected value: {value:.1f}")
            print(f"    EV advantage: {value - cash_out:.1f}")
            
        except Exception as e:
            print(f"  {description}: Error - {e}")
    
    # Test 8: Save Results
    print("\n8. Saving Results")
    print("-" * 40)
    
    final_results = {
        'qtable_analysis': analysis,
        'optimal_agent_results': optimal_results,
        'naive_comparison': {
            'naive_avg': naive_avg,
            'naive_std': naive_std,
            'naive_profitable': int(naive_profitable),
            'improvement': optimal_results['avg_reward'] - naive_avg,
            'improvement_percentage': ((optimal_results['avg_reward'] - naive_avg) / abs(naive_avg)) * 100
        }
    }
    
    with open('optimal_qtable_results.json', 'w') as f:
        # Convert numpy types for JSON serialization
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
    print("  - optimal_qtable.npz")
    print("  - optimal_qtable_results.json")
    
    print("\n" + "="*60)
    print("OPTIMAL Q-TABLE TEST COMPLETED")
    print("="*60)
    
    return final_results

if __name__ == "__main__":
    results = test_optimal_qtable()
    
    # Summary
    print(f"\nQUICK SUMMARY:")
    print(f"States with positive EV respins: {results['qtable_analysis']['percentage_states_with_positive_ev']:.1f}%")
    print(f"Optimal strategy average reward: {results['optimal_agent_results']['avg_reward']:.3f}")
    print(f"Improvement over naive: {results['naive_comparison']['improvement']:.3f} ({results['naive_comparison']['improvement_percentage']:.1f}%)")
    
    if results['optimal_agent_results']['profitable_rate'] > 0:
        print(f"Profitable episodes: {results['optimal_agent_results']['profitable_rate']:.1f}%")
    else:
        print("No profitable episodes found.")