import numpy as np
import random
from typing import Tuple, List
from tqdm import tqdm

class QLearningAgent:
    def __init__(self, 
                 state_space_size: int, 
                 action_space_size: int,
                 learning_rate: float = 0.1,
                 discount_factor: float = 0.95,
                 epsilon: float = 1.0,
                 epsilon_decay: float = 0.995,
                 epsilon_min: float = 0.01,
                 lr_decay: float = 0.9999,
                 lr_min: float = 0.01):
        
        self.state_space_size = state_space_size
        self.action_space_size = action_space_size
        self.learning_rate = learning_rate
        self.initial_learning_rate = learning_rate
        self.lr_decay = lr_decay
        self.lr_min = lr_min
        self.discount_factor = discount_factor
        self.epsilon = epsilon
        self.epsilon_decay = epsilon_decay
        self.epsilon_min = epsilon_min
        
        # Initialize Q-table with optimistic values to encourage exploration
        # Most actions will be slightly negative, but initialize optimistically
        self.q_table = np.random.normal(0.1, 0.05, (state_space_size, action_space_size))
        self.episode_rewards = []
        self.episode_count = 0
        
        # Performance tracking for adaptive learning
        self.recent_performance = []
        self.performance_window = 500
        self.best_recent_avg = float('-inf')
        
        # State importance tracking
        self.state_visit_counts = np.zeros(state_space_size)
        self.state_importance = np.ones(state_space_size)  # Higher values = more important states
        
        # Strategy consistency tracking
        self.state_action_preferences = np.zeros((state_space_size, action_space_size))
        self.strategy_consistency_bonus = 0.02
        
        # Decision temperature for progressive focusing
        self.initial_temperature = 2.0
        self.final_temperature = 0.5
        self.temperature_decay = 0.99995
    
    def choose_action(self, state_index: int, valid_actions: List[int] = None) -> int:
        if valid_actions is None:
            valid_actions = list(range(self.action_space_size))
        
        if random.random() < self.epsilon:
            # Random choice from valid actions only
            return random.choice(valid_actions)
        else:
            # Choose best action with more decisive selection
            q_values = self.q_table[state_index]
            valid_q_values = [q_values[action] for action in valid_actions]
            
            # Use adaptive temperature for progressive focusing
            current_temperature = max(self.final_temperature, 
                                    self.initial_temperature * (self.temperature_decay ** self.episode_count))
            valid_q_values = np.array(valid_q_values)
            
            # Normalize to prevent overflow
            max_q = np.max(valid_q_values)
            normalized_q = (valid_q_values - max_q) / current_temperature
            
            # Compute softmax safely
            exp_values = np.exp(normalized_q)
            probabilities = exp_values / np.sum(exp_values)
            
            # Handle any numerical issues
            if np.any(np.isnan(probabilities)) or np.sum(probabilities) == 0:
                # Fallback to greedy selection
                best_idx = np.argmax(valid_q_values)
                return valid_actions[best_idx]
            
            # Choose action based on softmax probabilities
            chosen_idx = np.random.choice(len(valid_actions), p=probabilities)
            return valid_actions[chosen_idx]
    
    def update_q_value(self, 
                      state_index: int, 
                      action: int, 
                      reward: float, 
                      next_state_index: int,
                      done: bool = False):
        
        # Track state visits and importance
        self.state_visit_counts[state_index] += 1
        self.state_action_preferences[state_index, action] += 1
        
        # Calculate importance-weighted learning rate
        importance_factor = self.state_importance[state_index]
        effective_lr = self.learning_rate * importance_factor
        
        current_q = self.q_table[state_index, action]
        max_next_q = np.max(self.q_table[next_state_index]) * (not done)
        
        td_error = reward + self.discount_factor * max_next_q - current_q
        new_q = current_q + effective_lr * td_error
        
        # Update state importance based on TD error magnitude
        # States with large TD errors are more important to learn from
        self.state_importance[state_index] = 0.99 * self.state_importance[state_index] + 0.01 * (1.0 + abs(td_error))
        
        self.q_table[state_index, action] = new_q
        
        # Encourage decisive Q-values: penalize states where all Q-values are too similar
        if self.episode_count > 10000:  # Only after some learning
            q_row = self.q_table[state_index]
            q_std = np.std(q_row)
            if q_std < 0.1:  # Q-values are too similar
                # Add small random perturbations to encourage differentiation
                self.q_table[state_index] += np.random.normal(0, 0.01, self.action_space_size)
    
    def decay_epsilon(self):
        """Conservative epsilon decay to maintain exploration."""
        if self.epsilon > self.epsilon_min:
            # Much more conservative decay - keep exploring longer
            if self.episode_count < 50000:
                # Very slow decay in first half of training
                decay_rate = 0.99992
            else:
                # Still slow decay in second half
                decay_rate = 0.99985
                
            self.epsilon *= decay_rate
    
    def decay_learning_rate(self):
        # Keep learning rate higher for longer to allow continued learning
        if self.learning_rate > self.lr_min and self.episode_count > 20000:
            self.learning_rate *= 0.99995  # Much slower decay
    
    def add_reward_shaping(self, state_index: int, action: int, game_state) -> float:
        """Add smart reward shaping based on game understanding."""
        bonus = 0.0
        
        # Get current symbol counts
        counts = game_state.get_current_counts()
        blank_count = counts[0]   # Blank
        coin_count = counts[1]    # Coin
        stack_count = counts[2]   # Coin Stack
        snake_count = counts[3]   # Snake
        net_count = counts[4]     # Net
        x2_count = counts[5]      # x2
        clover_count = counts[6]  # Clover
        crown_count = counts[7]   # Crown
        
        # Strong rewards for cashing out in good situations
        if action == 8:  # Cash out
            if clover_count > 0:
                bonus += 0.15 * clover_count  # Clovers are very valuable
            if snake_count > 0 and net_count == 0:
                bonus += 0.25 * snake_count   # Snakes without nets are great
            if crown_count >= 2:
                bonus += 0.2 * crown_count    # Multiple crowns are good
            if x2_count > 0 and (coin_count + stack_count) > 0:
                bonus += 0.1  # x2 with coins is valuable
                
        # Smart respin decisions
        else:
            respins_used = game_state.respins_used
            
            # Encourage respinning blanks early in the game
            if action == 0 and blank_count > 0 and respins_used < 3:
                bonus += 0.05
                
            # Encourage respinning nets when you have snakes
            if action == 4 and net_count > 0 and snake_count > 0:
                bonus += 0.15
                
            # Discourage risky respins when you already have good symbols
            if clover_count >= 2 and respins_used >= 3:
                bonus -= 0.1  # Don't risk good hands late
                
            # Encourage going for crown combinations
            if action == 7 and crown_count == 1 and respins_used < 4:
                bonus += 0.05  # Try for multiple crowns
                
            # Discourage excessive respinning
            if respins_used >= 4:
                bonus -= 0.05 * (respins_used - 3)
                
        # Strategy consistency bonus
        if self.state_visit_counts[state_index] > 5:
            # Reward consistent action choices in similar states
            total_visits = self.state_visit_counts[state_index]
            action_frequency = self.state_action_preferences[state_index, action] / total_visits
            
            # Bonus for taking the most preferred action in this state
            if action_frequency > 0.5:  # If this action is taken >50% of the time in this state
                bonus += self.strategy_consistency_bonus
                
        return bonus
    
    def train_episode(self, game, max_steps: int = 20) -> float:
        state = game.reset()
        state_index = game.get_current_abstract_state_index()  # Use abstract state
        total_reward = 0
        
        for step in range(max_steps):
            valid_actions = game.get_valid_actions()
            action = self.choose_action(state_index, valid_actions)
            reward, next_state, done = game.step(action)
            next_state_index = game.get_current_abstract_state_index()  # Use abstract state
            
            # Add reward shaping to help learning
            shaped_reward = reward + self.add_reward_shaping(state_index, action, game)
            self.update_q_value(state_index, action, shaped_reward, next_state_index, done)
            
            total_reward += reward
            state_index = next_state_index
            
            if done:
                break
        
        # Update performance tracking
        self.recent_performance.append(total_reward)
        if len(self.recent_performance) > self.performance_window:
            self.recent_performance.pop(0)
            
        # Track best recent performance for adaptive learning
        if len(self.recent_performance) >= 100:
            current_avg = np.mean(self.recent_performance[-100:])
            if current_avg > self.best_recent_avg:
                self.best_recent_avg = current_avg
                
            # Detect performance degradation and increase exploration
            if len(self.recent_performance) >= 500:
                recent_500 = np.mean(self.recent_performance[-500:])
                early_500 = np.mean(self.recent_performance[:500]) if len(self.recent_performance) >= 1000 else recent_500
                
                # If significantly worse than early performance, increase exploration
                if recent_500 < early_500 - 1.0 and self.epsilon < 0.1:
                    self.epsilon = min(0.2, self.epsilon * 2.0)  # Boost exploration
        
        self.decay_epsilon()
        self.decay_learning_rate()
        self.episode_rewards.append(total_reward)
        self.episode_count += 1
        
        return total_reward
    
    def train(self, game, num_episodes: int, max_steps_per_episode: int = 20, 
              verbose: bool = True, report_interval: int = 1000) -> List[float]:
        
        print(f"Starting training for {num_episodes} episodes...")
        print(f"Initial epsilon: {self.epsilon:.3f}")
        print(f"Learning rate: {self.learning_rate}")
        print(f"Discount factor: {self.discount_factor}")
        print("-" * 50)
        
        # Create progress bar
        pbar = tqdm(range(num_episodes), desc="Training", 
                   bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')
        
        for episode in pbar:
            episode_reward = self.train_episode(game, max_steps_per_episode)
            
            # Update progress bar description with current stats
            if (episode + 1) % 100 == 0:  # Update every 100 episodes to avoid too frequent updates
                recent_avg = np.mean(self.episode_rewards[-100:])
                pbar.set_postfix({
                    'avg_reward': f'{recent_avg:.2f}',
                    'epsilon': f'{self.epsilon:.3f}',
                    'lr': f'{self.learning_rate:.3f}',
                    'last_reward': f'{episode_reward:.2f}'
                })
            
            if verbose and (episode + 1) % report_interval == 0:
                avg_reward = np.mean(self.episode_rewards[-report_interval:])
                reward_std = np.std(self.episode_rewards[-report_interval:])
                best_reward = max(self.episode_rewards[-report_interval:])
                tqdm.write(f"\nEpisode {episode + 1}/{num_episodes}")
                tqdm.write(f"  Average reward (last {report_interval}): {avg_reward:.3f} ± {reward_std:.3f}")
                tqdm.write(f"  Best reward (last {report_interval}): {best_reward:.2f}")
                tqdm.write(f"  Current epsilon: {self.epsilon:.4f}")
                tqdm.write(f"  Current learning rate: {self.learning_rate:.4f}")
                tqdm.write(f"  Q-table sparsity: {np.count_nonzero(self.q_table) / self.q_table.size * 100:.1f}%")
        
        pbar.close()
        
        print(f"\nTraining completed!")
        print(f"Final epsilon: {self.epsilon:.4f}")
        print(f"Final learning rate: {self.learning_rate:.4f}")
        print(f"Average reward (last 1000 episodes): {np.mean(self.episode_rewards[-1000:]):.3f}")
        print(f"Best reward achieved: {max(self.episode_rewards):.2f}")
        print(f"Final Q-table sparsity: {np.count_nonzero(self.q_table) / self.q_table.size * 100:.1f}%")
        
        return self.episode_rewards
    
    def get_policy(self) -> np.ndarray:
        return np.argmax(self.q_table, axis=1)
    
    def evaluate_policy(self, game, num_episodes: int = 1000) -> Tuple[float, float]:
        old_epsilon = self.epsilon
        self.epsilon = 0.0  # No exploration during evaluation
        
        rewards = []
        
        # Create progress bar for evaluation
        pbar = tqdm(range(num_episodes), desc="Evaluating Policy", leave=False)
        
        for _ in pbar:
            state = game.reset()
            total_reward = 0
            
            for _ in range(20):  # Max steps per episode
                state_index = game.get_current_abstract_state_index()  # Use abstract state
                valid_actions = game.get_valid_actions()
                action = self.choose_action(state_index, valid_actions)
                reward, next_state, done = game.step(action)
                total_reward += reward
                
                if done:
                    break
            
            rewards.append(total_reward)
            
            # Update progress bar with running average
            if len(rewards) % 100 == 0:
                current_avg = np.mean(rewards)
                pbar.set_postfix({'avg_reward': f'{current_avg:.3f}'})
        
        pbar.close()
        
        self.epsilon = old_epsilon  # Restore original epsilon
        return np.mean(rewards), np.std(rewards)
    
    def get_state_values(self) -> np.ndarray:
        return np.max(self.q_table, axis=1)
    
    def analyze_strategy_focus(self) -> dict:
        """Analyze how focused/decisive the learned strategy is."""
        # Calculate Q-value standard deviation for each state
        q_stds = np.std(self.q_table, axis=1)
        
        # Calculate action entropy for each state (lower = more focused)
        action_entropies = []
        for state_idx in range(self.state_space_size):
            q_vals = self.q_table[state_idx]
            # Convert to probabilities using softmax
            exp_q = np.exp(q_vals - np.max(q_vals))  # Subtract max for numerical stability
            probs = exp_q / np.sum(exp_q)
            # Calculate entropy
            entropy = -np.sum(probs * np.log(probs + 1e-10))
            action_entropies.append(entropy)
        
        action_entropies = np.array(action_entropies)
        
        return {
            'mean_q_std': np.mean(q_stds),
            'mean_action_entropy': np.mean(action_entropies),
            'focused_states_pct': np.mean(q_stds > 0.5) * 100,  # States with clear preferences
            'decisive_states_pct': np.mean(action_entropies < 1.5) * 100  # Low entropy states
        }
    
    def save_model(self, filename: str):
        np.savez(filename, 
                 q_table=self.q_table,
                 episode_rewards=np.array(self.episode_rewards),
                 epsilon=self.epsilon,
                 episode_count=self.episode_count)
        print(f"Model saved to {filename}")
    
    def load_model(self, filename: str):
        data = np.load(filename)
        self.q_table = data['q_table']
        self.episode_rewards = data['episode_rewards'].tolist()
        self.epsilon = float(data['epsilon'])
        self.episode_count = int(data['episode_count'])
        print(f"Model loaded from {filename}")