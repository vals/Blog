import numpy as np
import random
from typing import Tuple, List
from tqdm import tqdm
from game_rules import SlotMachine
from expected_value_calculator import ExpectedValueCalculator

class HybridEVAgent:
    """Hybrid agent that uses expected value calculations to override default cash-out policy."""
    
    def __init__(self, game: SlotMachine, ev_threshold: float = 0.0):
        self.game = game
        self.ev_calculator = ExpectedValueCalculator(game)
        self.ev_threshold = ev_threshold  # Minimum EV advantage needed to override cash-out
        
        # Performance tracking
        self.episode_rewards = []
        self.episode_count = 0
        self.action_choices = {'cash_out': 0, 'ev_override': 0}
        self.ev_override_details = []
        
        # Strategy analysis
        self.state_action_log = []
        self.positive_ev_encounters = 0
        self.cash_out_decisions = 0
        
    def choose_action(self, state: Tuple[str, ...], respins_used: int, valid_actions: List[int] = None) -> int:
        """Choose action based on expected value calculations."""
        if valid_actions is None:
            valid_actions = self.game.get_valid_actions()
        
        # Calculate expected value for cashing out
        cash_out_value = self.ev_calculator._calculate_cash_out_value(state)
        
        # Find the best action and its expected value
        best_action, best_ev = self.ev_calculator.get_best_action(state, respins_used)
        
        # Check if we should override cash-out with a higher EV action
        if best_action != 8 and best_ev > cash_out_value + self.ev_threshold:
            # Override: Take the highest EV action
            self.action_choices['ev_override'] += 1
            self.positive_ev_encounters += 1
            
            # Log the override decision
            self.ev_override_details.append({
                'state': state,
                'respins_used': respins_used,
                'cash_out_value': cash_out_value,
                'chosen_action': best_action,
                'chosen_ev': best_ev,
                'ev_advantage': best_ev - cash_out_value
            })
            
            return best_action
        else:
            # Default: Cash out
            self.action_choices['cash_out'] += 1
            self.cash_out_decisions += 1
            return 8  # Cash out
    
    def simulate_episode(self, max_steps: int = 20) -> float:
        """Simulate one episode using the hybrid strategy."""
        state = self.game.reset()
        total_reward = 0
        
        for step in range(max_steps):
            valid_actions = self.game.get_valid_actions()
            action = self.choose_action(self.game.current_state, self.game.respins_used, valid_actions)
            
            # Log the state-action pair
            self.state_action_log.append({
                'state': self.game.current_state,
                'respins_used': self.game.respins_used,
                'action': action,
                'action_type': 'cash_out' if action == 8 else 'respin'
            })
            
            reward, next_state, done = self.game.step(action)
            total_reward += reward
            
            if done:
                break
        
        return total_reward
    
    def evaluate_strategy(self, num_episodes: int = 10000) -> dict:
        """Evaluate the hybrid strategy performance."""
        print(f"Evaluating Hybrid EV Strategy for {num_episodes} episodes...")
        
        rewards = []
        
        # Create progress bar
        pbar = tqdm(range(num_episodes), desc="Evaluating Strategy", 
                   bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')
        
        for episode in pbar:
            episode_reward = self.simulate_episode()
            rewards.append(episode_reward)
            self.episode_rewards.append(episode_reward)
            self.episode_count += 1
            
            # Update progress bar
            if (episode + 1) % 1000 == 0:
                recent_avg = np.mean(rewards[-1000:])
                pbar.set_postfix({
                    'avg_reward': f'{recent_avg:.2f}',
                    'ev_overrides': f'{self.action_choices["ev_override"]:.0f}',
                    'last_reward': f'{episode_reward:.2f}'
                })
        
        pbar.close()
        
        # Calculate performance metrics
        avg_reward = np.mean(rewards)
        std_reward = np.std(rewards)
        profitable_episodes = np.sum(np.array(rewards) > 0)
        best_reward = np.max(rewards)
        worst_reward = np.min(rewards)
        
        # Calculate action choice statistics
        total_decisions = sum(self.action_choices.values())
        cash_out_pct = (self.action_choices['cash_out'] / total_decisions) * 100
        ev_override_pct = (self.action_choices['ev_override'] / total_decisions) * 100
        
        results = {
            'avg_reward': avg_reward,
            'std_reward': std_reward,
            'profitable_episodes': profitable_episodes,
            'profitable_rate': (profitable_episodes / num_episodes) * 100,
            'best_reward': best_reward,
            'worst_reward': worst_reward,
            'cash_out_percentage': cash_out_pct,
            'ev_override_percentage': ev_override_pct,
            'total_decisions': total_decisions,
            'positive_ev_encounters': self.positive_ev_encounters,
            'ev_override_details': self.ev_override_details
        }
        
        return results
    
    def analyze_ev_overrides(self) -> dict:
        """Analyze the characteristics of EV override decisions."""
        if not self.ev_override_details:
            return {'message': 'No EV overrides occurred'}
        
        overrides = self.ev_override_details
        
        # Calculate statistics
        ev_advantages = [detail['ev_advantage'] for detail in overrides]
        cash_out_values = [detail['cash_out_value'] for detail in overrides]
        chosen_evs = [detail['chosen_ev'] for detail in overrides]
        
        # Analyze actions taken
        action_counts = {}
        for detail in overrides:
            action = detail['chosen_action']
            action_name = self.game.actions[action]['description']
            action_counts[action_name] = action_counts.get(action_name, 0) + 1
        
        # Analyze states where overrides occurred
        state_patterns = {}
        for detail in overrides:
            state = detail['state']
            # Create a pattern based on symbol counts
            counts = {}
            for symbol in self.game.symbols:
                counts[symbol] = state.count(symbol)
            
            # Find most common patterns
            pattern_key = f"Coins:{counts['Coin']}, Stacks:{counts['Coin Stack']}, Crowns:{counts['Crown']}, Clovers:{counts['Clover']}"
            state_patterns[pattern_key] = state_patterns.get(pattern_key, 0) + 1
        
        return {
            'total_overrides': len(overrides),
            'avg_ev_advantage': np.mean(ev_advantages),
            'max_ev_advantage': np.max(ev_advantages),
            'min_ev_advantage': np.min(ev_advantages),
            'avg_cash_out_value': np.mean(cash_out_values),
            'avg_chosen_ev': np.mean(chosen_evs),
            'action_distribution': action_counts,
            'common_state_patterns': dict(sorted(state_patterns.items(), key=lambda x: x[1], reverse=True)[:10])
        }
    
    def get_strategy_summary(self) -> dict:
        """Get a comprehensive summary of the strategy performance."""
        if not self.episode_rewards:
            return {'message': 'No episodes completed yet'}
        
        evaluation_results = {
            'total_episodes': self.episode_count,
            'avg_reward': np.mean(self.episode_rewards),
            'std_reward': np.std(self.episode_rewards),
            'profitable_episodes': np.sum(np.array(self.episode_rewards) > 0),
            'profitable_rate': (np.sum(np.array(self.episode_rewards) > 0) / len(self.episode_rewards)) * 100,
            'best_reward': np.max(self.episode_rewards),
            'worst_reward': np.min(self.episode_rewards)
        }
        
        decision_analysis = {
            'total_decisions': sum(self.action_choices.values()),
            'cash_out_decisions': self.action_choices['cash_out'],
            'ev_override_decisions': self.action_choices['ev_override'],
            'cash_out_percentage': (self.action_choices['cash_out'] / sum(self.action_choices.values())) * 100,
            'ev_override_percentage': (self.action_choices['ev_override'] / sum(self.action_choices.values())) * 100
        }
        
        return {
            'evaluation_results': evaluation_results,
            'decision_analysis': decision_analysis,
            'ev_override_analysis': self.analyze_ev_overrides()
        }
    
    def save_results(self, filename: str):
        """Save strategy results to file."""
        results = self.get_strategy_summary()
        
        # Also save raw data
        results['raw_data'] = {
            'episode_rewards': self.episode_rewards,
            'state_action_log': self.state_action_log,
            'ev_override_details': self.ev_override_details
        }
        
        np.savez(filename, **results)
        print(f"Hybrid EV strategy results saved to {filename}")
    
    def compare_with_naive_strategy(self, num_episodes: int = 10000) -> dict:
        """Compare hybrid strategy with naive cash-out-only strategy."""
        print("Comparing Hybrid EV Strategy with Naive Cash-Out Strategy...")
        
        # Run hybrid strategy
        hybrid_results = self.evaluate_strategy(num_episodes)
        
        # Run naive strategy (always cash out)
        naive_rewards = []
        for _ in tqdm(range(num_episodes), desc="Naive Strategy"):
            self.game.reset()
            reward, _, _ = self.game.step(8)  # Always cash out
            naive_rewards.append(reward)
        
        naive_avg = np.mean(naive_rewards)
        naive_std = np.std(naive_rewards)
        naive_profitable = np.sum(np.array(naive_rewards) > 0)
        
        comparison = {
            'hybrid_strategy': {
                'avg_reward': hybrid_results['avg_reward'],
                'std_reward': hybrid_results['std_reward'],
                'profitable_episodes': hybrid_results['profitable_episodes'],
                'profitable_rate': hybrid_results['profitable_rate']
            },
            'naive_strategy': {
                'avg_reward': naive_avg,
                'std_reward': naive_std,
                'profitable_episodes': naive_profitable,
                'profitable_rate': (naive_profitable / num_episodes) * 100
            },
            'performance_difference': {
                'reward_improvement': hybrid_results['avg_reward'] - naive_avg,
                'improvement_percentage': ((hybrid_results['avg_reward'] - naive_avg) / abs(naive_avg)) * 100,
                'profitable_rate_improvement': hybrid_results['profitable_rate'] - (naive_profitable / num_episodes) * 100
            }
        }
        
        return comparison