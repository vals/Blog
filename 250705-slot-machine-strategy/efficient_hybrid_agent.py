import numpy as np
import random
from typing import Tuple, List, Dict
from tqdm import tqdm
from game_rules import SlotMachine

class EfficientHybridAgent:
    """Efficient hybrid agent that uses pre-computed expected values for common scenarios."""
    
    def __init__(self, game: SlotMachine, ev_threshold: float = 0.0):
        self.game = game
        self.ev_threshold = ev_threshold
        self.symbol_probabilities = dict(zip(game.symbols, game.symbol_probabilities))
        
        # Performance tracking
        self.episode_rewards = []
        self.episode_count = 0
        self.action_choices = {'cash_out': 0, 'ev_override': 0}
        self.ev_override_details = []
        
        # Pre-compute common scenario EVs
        self.scenario_evs = self._precompute_scenario_evs()
        
    def _precompute_scenario_evs(self) -> Dict[str, float]:
        """Pre-compute expected values for common high-value scenarios."""
        evs = {}
        
        # Jackpot scenarios (3 Crowns -> try for 4th)
        crown_prob = self.symbol_probabilities['Crown']
        jackpot_ev = crown_prob * 100 - 1  # 8% chance of 100 Gold - 1 Gold cost
        evs['3_crowns_respin'] = jackpot_ev
        
        # 2 Crowns scenarios
        evs['2_crowns_respin'] = (crown_prob**2) * 100 - 2  # Need 2 more crowns
        
        # Coin completion scenarios
        coin_prob = self.symbol_probabilities['Coin']
        evs['2_coins_to_3'] = coin_prob * 3 - 1  # 30% chance of 3 Gold - 1 Gold cost
        evs['2_coins_to_4'] = (coin_prob**2) * 5 - 2  # Need 2 more coins for 5 Gold
        evs['3_coins_to_4'] = coin_prob * 5 - 1  # 30% chance of 5 Gold - 1 Gold cost
        
        # Coin Stack scenarios
        stack_prob = self.symbol_probabilities['Coin Stack']
        evs['2_stacks_to_3'] = stack_prob * 9 - 1  # 10% chance of 9 Gold - 1 Gold cost
        evs['2_stacks_to_4'] = (stack_prob**2) * 15 - 2  # Need 2 more stacks
        evs['3_stacks_to_4'] = stack_prob * 15 - 1  # 10% chance of 15 Gold - 1 Gold cost
        
        # Clover scenarios
        clover_prob = self.symbol_probabilities['Clover']
        evs['add_clover'] = clover_prob * 10 - 1  # 1% chance of 10 Gold - 1 Gold cost
        
        # Snake-Net scenarios
        snake_prob = self.symbol_probabilities['Snake']
        evs['add_snake_with_net'] = snake_prob * 3 - 1  # 10% chance of 3 Gold per net - 1 Gold cost
        
        return evs
    
    def _analyze_state_quickly(self, state: Tuple[str, ...]) -> Dict[str, int]:
        """Quickly analyze state for symbol counts."""
        counts = {}
        for symbol in self.game.symbols:
            counts[symbol] = state.count(symbol)
        return counts
    
    def _get_best_ev_action(self, state: Tuple[str, ...], respins_used: int) -> Tuple[int, float]:
        """Get best expected value action using pre-computed scenarios."""
        if respins_used >= self.game.max_respins:
            return 8, 0.0  # Must cash out
        
        counts = self._analyze_state_quickly(state)
        cash_out_value = self.game._calculate_payout(state, 1)
        
        best_action = 8
        best_ev = cash_out_value
        
        # Check high-value scenarios
        
        # 3 Crowns scenario - highest priority
        if counts['Crown'] == 3:
            respin_ev = cash_out_value + self.scenario_evs['3_crowns_respin']
            if respin_ev > best_ev:
                best_ev = respin_ev
                # Find action to respin non-Crown symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Crown':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # 2 Crowns scenario
        elif counts['Crown'] == 2:
            respin_ev = cash_out_value + self.scenario_evs['2_crowns_respin']
            if respin_ev > best_ev:
                best_ev = respin_ev
                # Find action to respin non-Crown symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Crown':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # 3 Coins scenario
        if counts['Coin'] == 3:
            respin_ev = cash_out_value + self.scenario_evs['3_coins_to_4']
            if respin_ev > best_ev:
                best_ev = respin_ev
                # Find action to respin non-Coin symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Coin':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # 2 Coins scenario
        elif counts['Coin'] == 2:
            respin_ev = cash_out_value + self.scenario_evs['2_coins_to_3']
            if respin_ev > best_ev:
                best_ev = respin_ev
                # Find action to respin non-Coin symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Coin':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # 3 Coin Stacks scenario
        if counts['Coin Stack'] == 3:
            respin_ev = cash_out_value + self.scenario_evs['3_stacks_to_4']
            if respin_ev > best_ev:
                best_ev = respin_ev
                # Find action to respin non-Coin Stack symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Coin Stack':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # 2 Coin Stacks scenario
        elif counts['Coin Stack'] == 2:
            respin_ev = cash_out_value + self.scenario_evs['2_stacks_to_3']
            if respin_ev > best_ev:
                best_ev = respin_ev
                # Find action to respin non-Coin Stack symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Coin Stack':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # Snake-Net interaction
        if counts['Net'] > 0 and counts['Snake'] < 4:
            # Try to add more snakes when we have nets
            add_snake_ev = cash_out_value + counts['Net'] * self.scenario_evs['add_snake_with_net']
            if add_snake_ev > best_ev:
                best_ev = add_snake_ev
                # Find action to respin non-Snake, non-Net symbol
                for i, symbol in enumerate(state):
                    if symbol not in ['Snake', 'Net']:
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # Clover scenarios - always try to add clovers if we have space
        if counts['Clover'] < 4:
            clover_ev = cash_out_value + self.scenario_evs['add_clover']
            if clover_ev > best_ev:
                best_ev = clover_ev
                # Find action to respin non-Clover symbol
                for i, symbol in enumerate(state):
                    if symbol != 'Clover':
                        symbol_action = self.game.symbols.index(symbol)
                        best_action = symbol_action
                        break
        
        # Remove snake if no net present
        if counts['Snake'] > 0 and counts['Net'] == 0:
            # Respin snakes to avoid zero payout
            snake_respin_ev = cash_out_value + 2  # Assume moderate benefit from removing snake
            if snake_respin_ev > best_ev:
                best_ev = snake_respin_ev
                best_action = 3  # Respin snake
        
        return best_action, best_ev
    
    def choose_action(self, state: Tuple[str, ...], respins_used: int, valid_actions: List[int] = None) -> int:
        """Choose action based on efficient expected value calculations."""
        if valid_actions is None:
            valid_actions = self.game.get_valid_actions()
        
        # Calculate cash-out value
        cash_out_value = self.game._calculate_payout(state, 1)
        
        # Get best EV action
        best_action, best_ev = self._get_best_ev_action(state, respins_used)
        
        # Check if we should override cash-out
        if best_action != 8 and best_ev > cash_out_value + self.ev_threshold:
            # Make sure the action is valid
            if best_action in valid_actions:
                self.action_choices['ev_override'] += 1
                self.ev_override_details.append({
                    'state': state,
                    'respins_used': respins_used,
                    'cash_out_value': cash_out_value,
                    'chosen_action': best_action,
                    'chosen_ev': best_ev,
                    'ev_advantage': best_ev - cash_out_value
                })
                return best_action
        
        # Default to cash out
        self.action_choices['cash_out'] += 1
        return 8
    
    def simulate_episode(self, max_steps: int = 20) -> float:
        """Simulate one episode using the efficient hybrid strategy."""
        state = self.game.reset()
        total_reward = 0
        
        for step in range(max_steps):
            valid_actions = self.game.get_valid_actions()
            action = self.choose_action(self.game.current_state, self.game.respins_used, valid_actions)
            
            reward, next_state, done = self.game.step(action)
            total_reward += reward
            
            if done:
                break
        
        return total_reward
    
    def evaluate_strategy(self, num_episodes: int = 10000) -> dict:
        """Evaluate the efficient hybrid strategy performance."""
        print(f"Evaluating Efficient Hybrid EV Strategy for {num_episodes} episodes...")
        
        rewards = []
        
        # Create progress bar
        pbar = tqdm(range(num_episodes), desc="Evaluating Strategy")
        
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
                    'cash_outs': f'{self.action_choices["cash_out"]:.0f}'
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
        cash_out_pct = (self.action_choices['cash_out'] / total_decisions) * 100 if total_decisions > 0 else 0
        ev_override_pct = (self.action_choices['ev_override'] / total_decisions) * 100 if total_decisions > 0 else 0
        
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
            'total_overrides': self.action_choices['ev_override'],
            'scenario_evs': self.scenario_evs
        }
        
        return results
    
    def compare_with_naive_strategy(self, num_episodes: int = 10000) -> dict:
        """Compare efficient hybrid strategy with naive cash-out-only strategy."""
        print(f"Comparing Efficient Hybrid EV Strategy with Naive Cash-Out Strategy...")
        
        # Run hybrid strategy
        hybrid_results = self.evaluate_strategy(num_episodes)
        
        # Run naive strategy (always cash out)
        naive_rewards = []
        print("Running Naive Strategy...")
        for _ in tqdm(range(num_episodes), desc="Naive Strategy"):
            self.game.reset()
            reward, _, _ = self.game.step(8)  # Always cash out
            naive_rewards.append(reward)
        
        naive_avg = np.mean(naive_rewards)
        naive_std = np.std(naive_rewards)
        naive_profitable = np.sum(np.array(naive_rewards) > 0)
        
        comparison = {
            'hybrid_strategy': hybrid_results,
            'naive_strategy': {
                'avg_reward': naive_avg,
                'std_reward': naive_std,
                'profitable_episodes': naive_profitable,
                'profitable_rate': (naive_profitable / num_episodes) * 100
            },
            'performance_difference': {
                'reward_improvement': hybrid_results['avg_reward'] - naive_avg,
                'improvement_percentage': ((hybrid_results['avg_reward'] - naive_avg) / abs(naive_avg)) * 100 if naive_avg != 0 else 0,
                'profitable_rate_improvement': hybrid_results['profitable_rate'] - (naive_profitable / num_episodes) * 100
            }
        }
        
        return comparison
    
    def get_override_analysis(self) -> dict:
        """Analyze the EV override decisions."""
        if not self.ev_override_details:
            return {'message': 'No EV overrides occurred'}
        
        overrides = self.ev_override_details
        ev_advantages = [detail['ev_advantage'] for detail in overrides]
        
        # Analyze actions taken
        action_counts = {}
        for detail in overrides:
            action = detail['chosen_action']
            action_name = self.game.actions.get(action, {}).get('description', f'Action {action}')
            action_counts[action_name] = action_counts.get(action_name, 0) + 1
        
        return {
            'total_overrides': len(overrides),
            'avg_ev_advantage': np.mean(ev_advantages),
            'max_ev_advantage': np.max(ev_advantages),
            'min_ev_advantage': np.min(ev_advantages),
            'action_distribution': action_counts
        }