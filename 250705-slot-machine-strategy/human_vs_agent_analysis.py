#!/usr/bin/env python3

"""
Human vs Q-Learning Agent Strategy Comparison Analysis

This script estimates a Q-table from empirical human data and compares it with 
the learned Q-learning agent strategy across multiple dimensions:
- Direct strategy comparison
- Behavioral pattern analysis  
- Value function insights
- Cognitive bias detection
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Tuple, List, Optional
from collections import defaultdict
import json
from game_rules import SlotMachine
from q_learning import QLearningAgent

class HumanQTableEstimator:
    """Estimates Q-table from empirical human decision data."""
    
    def __init__(self, alpha: float = 0.1, gamma: float = 0.99):
        self.alpha = alpha  # Learning rate for Q-table updates
        self.gamma = gamma  # Discount factor
        self.game = SlotMachine()
        self.q_table = {}  # Will store Q(s,a) values
        self.state_visits = defaultdict(int)
        self.state_action_visits = defaultdict(int)
        
    def state_to_abstract_index(self, symbol_counts: Dict[str, int], respins_used: int) -> int:
        """Convert symbol counts and respins to abstract state index (matching game abstraction)."""
        # Convert to tuple format expected by game
        counts_tuple = []
        for symbol in self.game.symbols:
            counts_tuple.append(symbol_counts.get(symbol, 0))
        counts_tuple.append(respins_used)
        
        return self.game.counts_to_index(tuple(counts_tuple))
    
    def _create_state_description(self, symbol_counts: Dict[str, int], respins_used: int) -> str:
        """Create human-readable description of game state."""
        # Create compact representation of non-zero symbols
        symbols = []
        for symbol, count in symbol_counts.items():
            if count > 0:
                if count == 1:
                    symbols.append(f"{symbol}")
                else:
                    symbols.append(f"{count}×{symbol}")
        
        if not symbols:
            symbol_desc = "Empty reels"
        else:
            symbol_desc = ", ".join(symbols)
        
        return f"{symbol_desc} (Respins: {respins_used})"
    
    def estimate_q_table_from_data(self, human_data: pd.DataFrame) -> Dict[int, Dict[int, float]]:
        """
        Estimate Q-table from human decision data using temporal difference learning.
        """
        print("🧠 Estimating Q-table from human decision data...")
        
        # Initialize Q-table and state descriptions
        self.state_descriptions = {}  # Store human-readable state descriptions
        
        for state in range(self.game.get_abstract_state_space_size()):
            self.q_table[state] = {}
            for action in range(self.game.get_action_space_size()):
                self.q_table[state][action] = 0.0
        
        # Group by episode to get sequential decisions
        episodes = human_data.groupby(['runId', 'episodeId'])
        
        for (run_id, episode_id), episode_data in episodes:
            episode_data = episode_data.sort_values('decisionId')
            
            for i in range(len(episode_data)):
                row = episode_data.iloc[i]
                
                # Get current state
                symbol_counts = {
                    'Blank': row['blank'], 'Coin': row['coin'], 
                    'Coin Stack': row['coinStack'], 'Snake': row['snake'],
                    'Net': row['net'], 'x2': row['x2'], 
                    'Clover': row['clover'], 'Crown': row['crown']
                }
                
                current_state = self.state_to_abstract_index(symbol_counts, row['respinsUsed'])
                action = int(row['action'])
                
                # Store human-readable state description
                if current_state not in self.state_descriptions:
                    self.state_descriptions[current_state] = self._create_state_description(symbol_counts, row['respinsUsed'])
                
                # Track visits
                self.state_visits[current_state] += 1
                self.state_action_visits[(current_state, action)] += 1
                
                # Get reward (calculate from balance change or payout)
                if i < len(episode_data) - 1:
                    next_row = episode_data.iloc[i + 1]
                    reward = next_row['balance'] - row['balance']
                    
                    # Get next state
                    next_symbol_counts = {
                        'Blank': next_row['blank'], 'Coin': next_row['coin'], 
                        'Coin Stack': next_row['coinStack'], 'Snake': next_row['snake'],
                        'Net': next_row['net'], 'x2': next_row['x2'], 
                        'Clover': next_row['clover'], 'Crown': next_row['crown']
                    }
                    next_state = self.state_to_abstract_index(next_symbol_counts, next_row['respinsUsed'])
                    
                    # Q-learning update: Q(s,a) += α[r + γ*max(Q(s',a')) - Q(s,a)]
                    max_next_q = max(self.q_table[next_state].values())
                    td_error = reward + self.gamma * max_next_q - self.q_table[current_state][action]
                    self.q_table[current_state][action] += self.alpha * td_error
                    
                else:
                    # Terminal state - reward is final payout
                    reward = row['currentPayout'] if row['action'] == 8 else 0
                    td_error = reward - self.q_table[current_state][action]
                    self.q_table[current_state][action] += self.alpha * td_error
        
        print(f"✅ Estimated Q-table from {len(episodes)} human episodes")
        print(f"   States visited: {len(self.state_visits)}")
        print(f"   State-action pairs: {len(self.state_action_visits)}")
        
        return self.q_table

class StrategyComparator:
    """Compares human and agent Q-learning strategies."""
    
    def __init__(self, human_q_table: Dict, agent_q_table: np.ndarray, human_data: pd.DataFrame, agent_data: pd.DataFrame, state_descriptions: Dict = None):
        self.human_q_table = human_q_table
        self.agent_q_table = agent_q_table
        self.human_data = human_data
        self.agent_data = agent_data
        self.game = SlotMachine()
        self.state_descriptions = state_descriptions or {}
        self.results = {}
        
    def compare_q_values(self) -> Dict:
        """Compare Q-values for identical states between human and agent."""
        print("📊 Comparing Q-values for identical states...")
        
        q_value_comparisons = []
        policy_agreements = []
        
        for state in range(self.agent_q_table.shape[0]):
            if state in self.human_q_table:
                human_q_values = [self.human_q_table[state][a] for a in range(9)]
                agent_q_values = self.agent_q_table[state].tolist()
                
                # Find best actions for each strategy
                human_best_action = np.argmax(human_q_values)
                agent_best_action = np.argmax(agent_q_values)
                
                policy_agreements.append(human_best_action == agent_best_action)
                
                # Calculate correlation between Q-values
                correlation = np.corrcoef(human_q_values, agent_q_values)[0, 1]
                if not np.isnan(correlation):
                    q_value_comparisons.append({
                        'state': state,
                        'human_best_action': human_best_action,
                        'agent_best_action': agent_best_action,
                        'agreement': human_best_action == agent_best_action,
                        'correlation': correlation,
                        'human_q_values': human_q_values,
                        'agent_q_values': agent_q_values
                    })
        
        policy_agreement_rate = np.mean(policy_agreements) * 100
        
        return {
            'policy_agreement_rate': policy_agreement_rate,
            'q_value_comparisons': q_value_comparisons,
            'states_compared': len(q_value_comparisons)
        }
    
    def analyze_risk_tolerance(self) -> Dict:
        """Analyze risk tolerance patterns in cash-out decisions."""
        print("🎲 Analyzing risk tolerance patterns...")
        
        def get_cash_out_rate_by_balance(data, label):
            cash_out_decisions = data[data['action'].isin([8])]  # Action 8 is cash out
            all_decisions = data.copy()
            
            # Group by balance ranges
            balance_ranges = [(0, 20), (21, 50), (51, 80), (81, 100), (101, float('inf'))]
            cash_out_rates = []
            
            for min_bal, max_bal in balance_ranges:
                range_decisions = all_decisions[
                    (all_decisions['balance'] >= min_bal) & 
                    (all_decisions['balance'] <= max_bal)
                ]
                range_cash_outs = cash_out_decisions[
                    (cash_out_decisions['balance'] >= min_bal) & 
                    (cash_out_decisions['balance'] <= max_bal)
                ]
                
                if len(range_decisions) > 0:
                    cash_out_rate = len(range_cash_outs) / len(range_decisions)
                    cash_out_rates.append({
                        'strategy': label,
                        'balance_range': f"{min_bal}-{max_bal if max_bal != float('inf') else '100+'}",
                        'cash_out_rate': cash_out_rate,
                        'total_decisions': len(range_decisions)
                    })
            
            return cash_out_rates
        
        human_risk = get_cash_out_rate_by_balance(self.human_data, 'Human')
        agent_risk = get_cash_out_rate_by_balance(self.agent_data, 'Agent')
        
        return {
            'human_risk_profile': human_risk,
            'agent_risk_profile': agent_risk
        }
    
    def analyze_special_symbol_handling(self) -> Dict:
        """Analyze how each strategy handles special symbols."""
        print("🐍 Analyzing special symbol handling...")
        
        def get_symbol_response_patterns(data, label):
            patterns = {}
            
            # Snake handling - look at decisions when snakes are present
            snake_situations = data[data['snake'] > 0]
            if len(snake_situations) > 0:
                snake_actions = snake_situations['action'].value_counts(normalize=True)
                patterns['snake_response'] = snake_actions.to_dict()
            
            # Net handling - especially with snakes present  
            net_snake_situations = data[(data['net'] > 0) & (data['snake'] > 0)]
            if len(net_snake_situations) > 0:
                net_snake_actions = net_snake_situations['action'].value_counts(normalize=True)
                patterns['net_snake_response'] = net_snake_actions.to_dict()
            
            # Multiplier (x2) handling
            multiplier_situations = data[data['x2'] > 0]
            if len(multiplier_situations) > 0:
                multiplier_actions = multiplier_situations['action'].value_counts(normalize=True)
                patterns['multiplier_response'] = multiplier_actions.to_dict()
            
            # Clover handling (high value symbol)
            clover_situations = data[data['clover'] > 0]
            if len(clover_situations) > 0:
                clover_actions = clover_situations['action'].value_counts(normalize=True)
                patterns['clover_response'] = clover_actions.to_dict()
            
            return patterns
        
        human_patterns = get_symbol_response_patterns(self.human_data, 'Human')
        agent_patterns = get_symbol_response_patterns(self.agent_data, 'Agent')
        
        return {
            'human_symbol_patterns': human_patterns,
            'agent_symbol_patterns': agent_patterns
        }
    
    def identify_divergent_states(self) -> Dict:
        """Identify states where strategies diverge most dramatically."""
        print("🔍 Identifying states with largest strategy divergences...")
        
        divergent_states = []
        
        for state in range(self.agent_q_table.shape[0]):
            if state in self.human_q_table:
                human_q_values = np.array([self.human_q_table[state][a] for a in range(9)])
                agent_q_values = self.agent_q_table[state]
                
                # Find best actions
                human_best = np.argmax(human_q_values)
                agent_best = np.argmax(agent_q_values)
                
                if human_best != agent_best:
                    # Calculate divergence magnitude
                    human_value_diff = human_q_values[human_best] - human_q_values[agent_best]
                    agent_value_diff = agent_q_values[agent_best] - agent_q_values[human_best]
                    
                    divergent_states.append({
                        'state': state,
                        'human_best_action': human_best,
                        'agent_best_action': agent_best,
                        'human_confidence': human_value_diff,
                        'agent_confidence': agent_value_diff,
                        'total_divergence': abs(human_value_diff) + abs(agent_value_diff)
                    })
        
        # Sort by total divergence
        divergent_states.sort(key=lambda x: x['total_divergence'], reverse=True)
        
        return {
            'most_divergent_states': divergent_states[:20],  # Top 20 most divergent
            'total_divergent_states': len(divergent_states)
        }
    
    def analyze_value_function_insights(self) -> Dict:
        """Compare which game states each strategy values most/least."""
        print("💎 Analyzing value function insights...")
        
        # Get max Q-values for each state (state values)
        human_state_values = []
        agent_state_values = []
        common_states = []
        
        for state in range(self.agent_q_table.shape[0]):
            if state in self.human_q_table:
                human_max_q = max(self.human_q_table[state].values())
                agent_max_q = np.max(self.agent_q_table[state])
                
                human_state_values.append(human_max_q)
                agent_state_values.append(agent_max_q)
                common_states.append(state)
        
        # Find states with highest/lowest values for each strategy
        human_values_array = np.array(human_state_values)
        agent_values_array = np.array(agent_state_values)
        
        # Top 10 most valued states for each strategy
        human_top_indices = np.argsort(human_values_array)[-10:]
        agent_top_indices = np.argsort(agent_values_array)[-10:]
        
        human_top_states = [common_states[i] for i in human_top_indices]
        agent_top_states = [common_states[i] for i in agent_top_indices]
        
        # Bottom 10 least valued states
        human_bottom_indices = np.argsort(human_values_array)[:10]
        agent_bottom_indices = np.argsort(agent_values_array)[:10]
        
        human_bottom_states = [common_states[i] for i in human_bottom_indices]
        agent_bottom_states = [common_states[i] for i in agent_bottom_indices]
        
        return {
            'human_most_valued_states': list(zip(human_top_states, human_values_array[human_top_indices])),
            'agent_most_valued_states': list(zip(agent_top_states, agent_values_array[agent_top_indices])),
            'human_least_valued_states': list(zip(human_bottom_states, human_values_array[human_bottom_indices])),
            'agent_least_valued_states': list(zip(agent_bottom_states, agent_values_array[agent_bottom_indices])),
            'value_correlation': np.corrcoef(human_values_array, agent_values_array)[0, 1]
        }
    
    def detect_cognitive_biases(self) -> Dict:
        """Detect potential cognitive biases in human decision patterns."""
        print("🧠 Detecting cognitive biases in human decisions...")
        
        biases = {}
        
        # Loss aversion: tendency to make riskier decisions when losing
        losing_decisions = self.human_data[self.human_data['balance'] < 100]
        winning_decisions = self.human_data[self.human_data['balance'] >= 100]
        
        if len(losing_decisions) > 0 and len(winning_decisions) > 0:
            losing_cash_out_rate = len(losing_decisions[losing_decisions['action'] == 8]) / len(losing_decisions)
            winning_cash_out_rate = len(winning_decisions[winning_decisions['action'] == 8]) / len(winning_decisions)
            
            biases['loss_aversion_indicator'] = {
                'losing_cash_out_rate': losing_cash_out_rate,
                'winning_cash_out_rate': winning_cash_out_rate,
                'difference': winning_cash_out_rate - losing_cash_out_rate
            }
        
        # Gambler's fallacy: changing behavior after wins/losses
        human_episodes = self.human_data.groupby(['runId', 'episodeId'])
        consecutive_patterns = []
        
        for (run_id, episode_id), episode_data in human_episodes:
            episode_data = episode_data.sort_values('decisionId')
            if len(episode_data) >= 2:
                # Look at decision patterns within episodes
                actions = episode_data['action'].values
                for i in range(len(actions) - 1):
                    consecutive_patterns.append((actions[i], actions[i+1]))
        
        if consecutive_patterns:
            # Simple pattern analysis
            pattern_counts = {}
            for pattern in consecutive_patterns:
                pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1
            
            biases['decision_sequence_patterns'] = dict(sorted(pattern_counts.items(), key=lambda x: x[1], reverse=True)[:10])
        
        # Overconfidence: continuing to play with low balances
        low_balance_decisions = self.human_data[self.human_data['balance'] <= 20]
        if len(low_balance_decisions) > 0:
            continue_rate_low_balance = 1 - (len(low_balance_decisions[low_balance_decisions['action'] == 8]) / len(low_balance_decisions))
            biases['overconfidence_indicator'] = {
                'continue_rate_at_low_balance': continue_rate_low_balance,
                'decisions_at_low_balance': len(low_balance_decisions)
            }
        
        return biases
    
    def generate_comparison_report(self) -> Dict:
        """Generate comprehensive comparison report."""
        print("\n🎯 Generating comprehensive strategy comparison report...")
        
        # Run all analyses
        q_value_comparison = self.compare_q_values()
        risk_analysis = self.analyze_risk_tolerance()
        symbol_analysis = self.analyze_special_symbol_handling()
        divergent_states = self.identify_divergent_states()
        value_insights = self.analyze_value_function_insights()
        bias_detection = self.detect_cognitive_biases()
        
        self.results = {
            'summary': {
                'total_states_compared': q_value_comparison['states_compared'],
                'policy_agreement_rate': q_value_comparison['policy_agreement_rate'],
                'value_correlation': value_insights['value_correlation'],
                'divergent_states_count': divergent_states['total_divergent_states']
            },
            'q_value_comparison': q_value_comparison,
            'risk_tolerance': risk_analysis,
            'special_symbol_handling': symbol_analysis,
            'divergent_states': divergent_states,
            'value_function_insights': value_insights,
            'cognitive_biases': bias_detection,
            'state_descriptions': self.state_descriptions
        }
        
        return self.results
    
    def save_results(self, filename: str = "human_vs_agent_analysis.json"):
        """Save analysis results to JSON file."""
        # Convert numpy types to native Python types for JSON serialization
        def convert_numpy(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.bool_):
                return bool(obj)
            return obj
        
        def clean_for_json(data):
            if isinstance(data, dict):
                cleaned_dict = {}
                for k, v in data.items():
                    # Convert tuple keys to strings
                    if isinstance(k, tuple):
                        key = str(k)
                    else:
                        key = convert_numpy(k)
                    cleaned_dict[key] = clean_for_json(v)
                return cleaned_dict
            elif isinstance(data, list):
                return [clean_for_json(item) for item in data]
            else:
                return convert_numpy(data)
        
        clean_results = clean_for_json(self.results)
        
        with open(filename, 'w') as f:
            json.dump(clean_results, f, indent=2)
        
        print(f"📊 Analysis results saved to {filename}")

def main():
    """Main analysis function."""
    print("🎰 Human vs Q-Learning Agent Strategy Comparison")
    print("=" * 60)
    
    # Load data
    print("📚 Loading datasets...")
    human_data = pd.read_csv("slot_machine_data.csv")
    agent_data = pd.read_csv("simulated_slot_machine_data.csv")
    
    print(f"Human data: {len(human_data)} decisions from {human_data['episodeId'].nunique()} episodes")
    print(f"Agent data: {len(agent_data)} decisions from {agent_data['episodeId'].nunique()} episodes")
    
    # Estimate human Q-table
    estimator = HumanQTableEstimator()
    human_q_table = estimator.estimate_q_table_from_data(human_data)
    
    # Load agent Q-table
    print("🤖 Loading agent Q-table...")
    agent = QLearningAgent(
        state_space_size=SlotMachine().get_abstract_state_space_size(),
        action_space_size=SlotMachine().get_action_space_size(),
        learning_rate=0.25,
        discount_factor=0.99,
        epsilon=0.0
    )
    agent.load_model("slot_machine_model.npz")
    agent_q_table = agent.q_table
    
    # Create comprehensive state descriptions (including agent-only states)
    all_state_descriptions = estimator.state_descriptions.copy()
    
    # Add descriptions for agent-only states by looking at agent data
    print("🔍 Creating descriptions for agent-only states...")
    agent_states_seen = set()
    agent_episodes = agent_data.groupby(['runId', 'episodeId'])
    
    for (run_id, episode_id), episode_data in agent_episodes:
        episode_data = episode_data.sort_values('decisionId')
        for _, row in episode_data.iterrows():
            symbol_counts = {
                'Blank': row['blank'], 'Coin': row['coin'], 
                'Coin Stack': row['coinStack'], 'Snake': row['snake'],
                'Net': row['net'], 'x2': row['x2'], 
                'Clover': row['clover'], 'Crown': row['crown']
            }
            state_index = estimator.state_to_abstract_index(symbol_counts, row['respinsUsed'])
            if state_index not in all_state_descriptions:
                all_state_descriptions[state_index] = estimator._create_state_description(symbol_counts, row['respinsUsed'])
            agent_states_seen.add(state_index)
    
    print(f"   Added descriptions for {len(all_state_descriptions) - len(estimator.state_descriptions)} additional states")
    
    # Try to generate descriptions for remaining high-value states by brute force
    print("🧮 Attempting to decode remaining high-value states...")
    
    def generate_possible_states():
        """Generate possible symbol combinations for common game situations."""
        possible_states = []
        symbols = ['Blank', 'Coin', 'Coin Stack', 'Snake', 'Net', 'x2', 'Clover', 'Crown']
        
        # Generate some high-value combinations that might be missing
        high_value_combinations = [
            # Jackpot scenarios
            ({'Crown': 4}, 0), ({'Crown': 4}, 1), ({'Crown': 4}, 2), ({'Crown': 4}, 4),
            # High coin scenarios  
            ({'Coin': 4}, 0), ({'Coin': 4}, 1), ({'Coin': 4}, 2),
            ({'Coin': 3, 'x2': 1}, 0), ({'Coin': 3, 'x2': 1}, 1),
            # Coin stack scenarios
            ({'Coin Stack': 4}, 0), ({'Coin Stack': 4}, 1), ({'Coin Stack': 4}, 2),
            ({'Coin Stack': 3, 'x2': 1}, 0), ({'Coin Stack': 3, 'x2': 1}, 1),
            # Clover scenarios
            ({'Clover': 1, 'Coin': 3}, 0), ({'Clover': 2, 'Coin': 2}, 0), ({'Clover': 3, 'Coin': 1}, 0),
            ({'Clover': 1, 'x2': 1, 'Coin': 2}, 0), ({'Clover': 1, 'Crown': 3}, 0),
            # Mixed high-value scenarios
            ({'Coin Stack': 2, 'Coin': 2}, 0), ({'Coin Stack': 2, 'Crown': 2}, 0),
            ({'Crown': 3, 'Coin': 1}, 0), ({'Crown': 3, 'x2': 1}, 0),
        ]
        
        for symbol_counts, respins in high_value_combinations:
            # Fill remaining slots with blanks
            total_symbols = sum(symbol_counts.values())
            if total_symbols <= 4:
                full_counts = {symbol: 0 for symbol in symbols}
                full_counts.update(symbol_counts)
                full_counts['Blank'] = 4 - total_symbols
                possible_states.append((full_counts, respins))
        
        return possible_states
    
    # Generate descriptions for missing states
    possible_states = generate_possible_states()
    for symbol_counts, respins in possible_states:
        state_index = estimator.state_to_abstract_index(symbol_counts, respins)
        if state_index not in all_state_descriptions:
            all_state_descriptions[state_index] = estimator._create_state_description(symbol_counts, respins)
    
    print(f"   Generated descriptions for {len(all_state_descriptions) - len(estimator.state_descriptions) - 23} additional theoretical states")
    
    # Run comparison analysis
    comparator = StrategyComparator(human_q_table, agent_q_table, human_data, agent_data, all_state_descriptions)
    results = comparator.generate_comparison_report()
    
    # Save results
    comparator.save_results()
    
    # Print summary
    print("\n🎯 ANALYSIS SUMMARY")
    print("=" * 40)
    print(f"States compared: {results['summary']['total_states_compared']}")
    print(f"Policy agreement rate: {results['summary']['policy_agreement_rate']:.1f}%")
    print(f"Value correlation: {results['summary']['value_correlation']:.3f}")
    print(f"Divergent states: {results['summary']['divergent_states_count']}")
    
    print("\n✅ Analysis complete! Check 'human_vs_agent_analysis.json' for detailed results.")

if __name__ == "__main__":
    main()