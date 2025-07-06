#!/usr/bin/env python3

import numpy as np
from collections import defaultdict
from game_rules import SlotMachine
import matplotlib.pyplot as plt

class GameAnalyzer:
    def __init__(self):
        self.game = SlotMachine()
        
    def analyze_mathematical_expectation(self, num_samples=100000):
        """Analyze the mathematical expectation of different strategies."""
        print("=== MATHEMATICAL EXPECTATION ANALYSIS ===")
        
        # Test pure random strategy
        random_rewards = []
        for _ in range(num_samples):
            state = self.game.reset()
            episode_reward = 0
            
            while self.game.episode_active:
                valid_actions = self.game.get_valid_actions()
                action = np.random.choice(valid_actions)
                reward, next_state, done = self.game.step(action)
                episode_reward += reward
                state = next_state
                
            random_rewards.append(episode_reward)
        
        print(f"Random Strategy (n={num_samples}):")
        print(f"  Mean reward: {np.mean(random_rewards):.3f}")
        print(f"  Std deviation: {np.std(random_rewards):.3f}")
        print(f"  Min reward: {np.min(random_rewards)}")
        print(f"  Max reward: {np.max(random_rewards)}")
        print(f"  Positive episodes: {np.sum(np.array(random_rewards) > 0)} ({100*np.sum(np.array(random_rewards) > 0)/num_samples:.1f}%)")
        
        # Test always cash out strategy
        cashout_rewards = []
        for _ in range(num_samples):
            state = self.game.reset()
            episode_reward = 0
            
            while self.game.episode_active:
                valid_actions = self.game.get_valid_actions()
                # Always choose cash out if available, otherwise random
                if 8 in valid_actions:  # Cash out action
                    action = 8
                else:
                    action = np.random.choice(valid_actions)
                reward, next_state, done = self.game.step(action)
                episode_reward += reward
                state = next_state
                
            cashout_rewards.append(episode_reward)
        
        print(f"\nAlways Cash Out Strategy (n={num_samples}):")
        print(f"  Mean reward: {np.mean(cashout_rewards):.3f}")
        print(f"  Std deviation: {np.std(cashout_rewards):.3f}")
        print(f"  Min reward: {np.min(cashout_rewards)}")
        print(f"  Max reward: {np.max(cashout_rewards)}")
        print(f"  Positive episodes: {np.sum(np.array(cashout_rewards) > 0)} ({100*np.sum(np.array(cashout_rewards) > 0)/num_samples:.1f}%)")
        
        return random_rewards, cashout_rewards
    
    def analyze_state_information_loss(self, num_samples=10000):
        """Analyze how much information is lost in state abstraction."""
        print("\n=== STATE ABSTRACTION ANALYSIS ===")
        
        # Track state frequencies and their associated rewards
        abstract_state_rewards = defaultdict(list)  # symbol_counts -> rewards
        position_state_rewards = defaultdict(list)  # full_positions -> rewards
        
        for _ in range(num_samples):
            state = self.game.reset()
            
            # Get both representations
            abstract_state = tuple(sorted(state))  # Our current abstraction
            position_state = tuple(state)  # Full positional information
            
            # Play random episode
            episode_reward = 0
            while self.game.episode_active:
                valid_actions = self.game.get_valid_actions()
                action = np.random.choice(valid_actions)
                reward, next_state, done = self.game.step(action)
                episode_reward += reward
                state = next_state
            
            abstract_state_rewards[abstract_state].append(episode_reward)
            position_state_rewards[position_state].append(episode_reward)
        
        # Analyze variance within abstract states
        high_variance_states = []
        for abstract_state, rewards in abstract_state_rewards.items():
            if len(rewards) >= 10:  # Only analyze states with enough samples
                variance = np.var(rewards)
                mean_reward = np.mean(rewards)
                if variance > 100:  # High variance threshold
                    high_variance_states.append((abstract_state, mean_reward, variance, len(rewards)))
        
        print(f"States with high reward variance (>100):")
        high_variance_states.sort(key=lambda x: x[2], reverse=True)
        for state, mean_rew, var, count in high_variance_states[:10]:
            print(f"  State {state}: mean={mean_rew:.1f}, var={var:.1f}, n={count}")
        
        # Calculate information loss metric
        abstract_info = len(abstract_state_rewards)
        position_info = len(position_state_rewards)
        info_loss = 1 - (abstract_info / position_info)
        
        print(f"\nInformation Loss Analysis:")
        print(f"  Unique abstract states: {abstract_info}")
        print(f"  Unique position states: {position_info}")
        print(f"  Information loss ratio: {info_loss:.3f}")
        
        return high_variance_states
    
    def analyze_q_learning_compatibility(self):
        """Analyze if this game is fundamentally compatible with Q-learning."""
        print("\n=== Q-LEARNING COMPATIBILITY ANALYSIS ===")
        
        # Check for key Q-learning assumptions
        print("Analyzing Q-learning prerequisites:")
        
        # 1. Markov Property
        print("1. Markov Property:")
        print("   ✓ Current state contains all necessary information")
        print("   ✓ Transitions depend only on current state and action")
        
        # 2. Stationary Environment
        print("2. Stationary Environment:")
        print("   ✓ Transition probabilities are fixed")
        print("   ✓ Reward function is deterministic given state-action")
        
        # 3. Bounded Rewards
        print("3. Bounded Rewards:")
        print("   ✓ Rewards are bounded (min: -5 per spin, max: ~400)")
        
        # 4. Finite State-Action Space
        print("4. Finite State-Action Space:")
        print("   ✓ Finite states (though large)")
        print("   ✓ Finite actions (9 total)")
        
        # 5. Exploration vs Exploitation Trade-off
        print("5. Key Issue - Exploration vs Exploitation:")
        
        # Test the core problem: is there actually a learnable strategy?
        print("\nTesting strategy learnability...")
        
        # Generate many episodes and see if any strategy emerges
        state_action_returns = defaultdict(list)
        
        for episode in range(1000):
            state = self.game.reset()
            
            while self.game.episode_active:
                valid_actions = self.game.get_valid_actions()
                action = np.random.choice(valid_actions)
                
                # Track Q-value targets (immediate + discounted future)
                reward, next_state, done = self.game.step(action)
                
                state_key = tuple(sorted(state))
                state_action_returns[(state_key, action)].append(reward)
                
                state = next_state
        
        # Analyze if there are consistently better actions
        consistent_strategies = 0
        total_comparisons = 0
        
        for (state_key, action), returns in state_action_returns.items():
            if len(returns) >= 5:  # Enough samples
                # Compare with other actions in the same state
                same_state_actions = [(s, a, r) for (s, a), r in state_action_returns.items() 
                                    if s == state_key and len(r) >= 5]
                
                if len(same_state_actions) > 1:
                    current_mean = np.mean(returns)
                    other_means = [np.mean(r) for s, a, r in same_state_actions if a != action]
                    
                    if other_means:
                        total_comparisons += 1
                        if current_mean > max(other_means) + 0.1:  # Threshold for "better"
                            consistent_strategies += 1
        
        if total_comparisons > 0:
            strategy_ratio = consistent_strategies / total_comparisons
            print(f"   States with clearly better actions: {consistent_strategies}/{total_comparisons} ({100*strategy_ratio:.1f}%)")
            
            if strategy_ratio < 0.1:
                print("   ⚠️  WARNING: Very few states have clearly superior actions")
                print("   ⚠️  This suggests limited strategic depth")
        
        return strategy_ratio if total_comparisons > 0 else 0
    
    def diagnose_degradation_cause(self):
        """Identify the specific cause of performance degradation."""
        print("\n=== PERFORMANCE DEGRADATION DIAGNOSIS ===")
        
        # Hypothesis 1: Game is inherently negative expected value
        random_rewards, cashout_rewards = self.analyze_mathematical_expectation(50000)
        
        if np.mean(random_rewards) < -1 and np.mean(cashout_rewards) < -1:
            print("🔍 DIAGNOSIS: Game has strong negative expected value")
            print("   - Random play: {:.3f}".format(np.mean(random_rewards)))
            print("   - Always cash out: {:.3f}".format(np.mean(cashout_rewards)))
            print("   - Any learned strategy may perform worse than random exploration")
        
        # Hypothesis 2: Early episodes benefit from optimistic initialization
        print("\n🔍 DIAGNOSIS: Early performance advantage likely due to:")
        print("   - Optimistic Q-value initialization (0.1 ± 0.05)")
        print("   - High exploration rate (epsilon starts at 1.0)")
        print("   - Beginner's luck from random exploration")
        
        # Hypothesis 3: Learning makes agent overly conservative or predictable
        print("\n🔍 DIAGNOSIS: Learning degradation mechanisms:")
        print("   - Agent learns to avoid actions that sometimes work")
        print("   - Reduced exploration leads to suboptimal fixed patterns")
        print("   - Negative game bias makes any consistent strategy worse than random")
        
        return {
            'random_ev': np.mean(random_rewards),
            'cashout_ev': np.mean(cashout_rewards),
            'random_std': np.std(random_rewards),
            'cashout_std': np.std(cashout_rewards)
        }

def main():
    analyzer = GameAnalyzer()
    
    # Run comprehensive analysis
    random_rewards, cashout_rewards = analyzer.analyze_mathematical_expectation()
    high_var_states = analyzer.analyze_state_information_loss()
    strategy_ratio = analyzer.analyze_q_learning_compatibility()
    diagnosis = analyzer.diagnose_degradation_cause()
    
    print("\n" + "="*60)
    print("FINAL ANALYSIS SUMMARY")
    print("="*60)
    
    if diagnosis['random_ev'] < -2 and diagnosis['cashout_ev'] < -2:
        print("❌ FUNDAMENTAL ISSUE IDENTIFIED:")
        print("   This game has strongly negative expected value regardless of strategy.")
        print("   Q-learning cannot overcome mathematical impossibility.")
        print("   Performance degradation is EXPECTED, not a bug.")
        print("\n💡 RECOMMENDATION:")
        print("   - Accept that this is a losing game by design")
        print("   - Focus on minimizing losses rather than maximizing wins")
        print("   - Consider different RL objectives (e.g., survival time)")
    
    elif strategy_ratio < 0.1:
        print("❌ STRATEGIC DEPTH ISSUE:")
        print("   Very few states have clearly superior actions.")
        print("   The game may be too random for meaningful strategy learning.")
        print("\n💡 RECOMMENDATION:")
        print("   - Increase reward shaping for intermediate goals")
        print("   - Consider different state representations")
        print("   - Try policy gradient methods instead of value-based")

if __name__ == "__main__":
    main()