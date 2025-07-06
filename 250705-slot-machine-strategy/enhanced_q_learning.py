import numpy as np
import random
from typing import Tuple, List
from tqdm import tqdm

class EnhancedQLearningAgent:
    def __init__(self, 
                 state_space_size: int, 
                 action_space_size: int,
                 learning_rate: float = 0.25,
                 discount_factor: float = 0.99,
                 epsilon: float = 1.0,
                 epsilon_decay: float = 0.9995,
                 epsilon_min: float = 0.01,
                 lr_decay: float = 0.9999,
                 lr_min: float = 0.01,
                 use_double_q: bool = True,
                 adaptive_epsilon: bool = True):
        
        self.state_space_size = state_space_size
        self.action_space_size = action_space_size
        self.learning_rate = learning_rate
        self.initial_learning_rate = learning_rate
        self.lr_decay = lr_decay
        self.lr_min = lr_min
        self.discount_factor = discount_factor
        self.epsilon = epsilon
        self.initial_epsilon = epsilon
        self.epsilon_decay = epsilon_decay
        self.epsilon_min = epsilon_min
        self.use_double_q = use_double_q
        self.adaptive_epsilon = adaptive_epsilon
        
        # Double Q-Learning: Two Q-tables
        if use_double_q:
            self.q_table_a = np.zeros((state_space_size, action_space_size))
            self.q_table_b = np.zeros((state_space_size, action_space_size))
        else:
            self.q_table = np.zeros((state_space_size, action_space_size))
        
        # Tracking
        self.episode_rewards = []
        self.episode_count = 0
        self.recent_performance = []
        self.performance_window = 1000
        
        # Enhanced tracking for adaptive epsilon
        self.q_value_changes = []
        self.exploration_bonus = 0.0
    
    def choose_action(self, state_index: int, valid_actions: List[int] = None) -> int:
        if valid_actions is None:
            valid_actions = list(range(self.action_space_size))
        
        if random.random() < self.epsilon:
            # Random choice from valid actions only
            return random.choice(valid_actions)
        else:
            # Choose best action from valid actions only
            if self.use_double_q:
                # Use average of both Q-tables for action selection
                combined_q = (self.q_table_a[state_index] + self.q_table_b[state_index]) / 2
            else:
                combined_q = self.q_table[state_index]
            
            masked_q_values = {action: combined_q[action] for action in valid_actions}
            return max(masked_q_values, key=masked_q_values.get)
    
    def update_q_value(self, 
                      state_index: int, 
                      action: int, 
                      reward: float, 
                      next_state_index: int):
        
        if self.use_double_q:
            # Double Q-Learning update
            if random.random() < 0.5:
                # Update Q_A, use Q_B for target
                old_q = self.q_table_a[state_index, action]
                best_next_action = np.argmax(self.q_table_a[next_state_index])
                max_next_q = self.q_table_b[next_state_index, best_next_action]
                
                new_q = old_q + self.learning_rate * (
                    reward + self.discount_factor * max_next_q - old_q
                )
                self.q_table_a[state_index, action] = new_q
            else:
                # Update Q_B, use Q_A for target
                old_q = self.q_table_b[state_index, action]
                best_next_action = np.argmax(self.q_table_b[next_state_index])
                max_next_q = self.q_table_a[next_state_index, best_next_action]
                
                new_q = old_q + self.learning_rate * (
                    reward + self.discount_factor * max_next_q - old_q
                )
                self.q_table_b[state_index, action] = new_q
        else:
            # Standard Q-Learning update
            old_q = self.q_table[state_index, action]
            max_next_q = np.max(self.q_table[next_state_index])
            
            new_q = old_q + self.learning_rate * (
                reward + self.discount_factor * max_next_q - old_q
            )
            self.q_table[state_index, action] = new_q
        
        # Track Q-value changes for adaptive learning
        q_change = abs(new_q - old_q)
        self.q_value_changes.append(q_change)
        if len(self.q_value_changes) > 1000:
            self.q_value_changes.pop(0)
    
    def decay_epsilon(self):
        if self.adaptive_epsilon and len(self.recent_performance) > 100:
            # Adaptive epsilon based on learning progress
            recent_improvement = self._calculate_improvement_rate()
            
            if recent_improvement > 0.1:  # Good learning progress
                # Decay faster when learning well
                decay_rate = self.epsilon_decay * 0.99
            elif recent_improvement < -0.05:  # Poor or negative progress
                # Decay slower or even increase epsilon
                decay_rate = self.epsilon_decay * 1.001
            else:
                # Normal decay
                decay_rate = self.epsilon_decay
            
            if self.epsilon > self.epsilon_min:
                self.epsilon *= decay_rate
        else:
            # Standard epsilon decay
            if self.epsilon > self.epsilon_min:
                self.epsilon *= self.epsilon_decay
    
    def _calculate_improvement_rate(self) -> float:
        """Calculate recent performance improvement rate."""
        if len(self.recent_performance) < 200:
            return 0.0
        
        # Compare recent vs older performance
        recent_avg = np.mean(self.recent_performance[-100:])
        older_avg = np.mean(self.recent_performance[-200:-100])
        
        if older_avg == 0:
            return 0.0
        
        return (recent_avg - older_avg) / abs(older_avg)
    
    def decay_learning_rate(self):
        if self.learning_rate > self.lr_min:
            self.learning_rate *= self.lr_decay
    
    def add_reward_shaping(self, state_index: int, action: int, game_state) -> float:
        """Add intermediate rewards for good decision making."""
        bonus = 0.0
        
        # Get symbol counts from game state
        counts = game_state.get_current_counts()
        snake_count = counts[3]  # Snake is index 3
        net_count = counts[4]    # Net is index 4
        clover_count = counts[6] # Clover is index 6
        
        # Reward smart decisions
        if action == 8:  # Cash out
            # Bonus for cashing out with good symbols
            if clover_count > 0:
                bonus += 0.1 * clover_count  # Reward keeping clovers
            if snake_count > 0 and net_count == 0:
                bonus += 0.2  # Smart to cash out with unprotected snake
        
        elif action == 3:  # Respin snake
            if net_count == 0:
                bonus += 0.1  # Good to remove dangerous snakes
        
        elif action == 0:  # Respin blank
            if clover_count > 0 or (snake_count > 0 and net_count > 0):
                bonus += 0.05  # Reasonable to try improving blanks
        
        return bonus
    
    def train_episode(self, game, max_steps: int = 20) -> float:
        state = game.reset()
        state_index = game.get_current_abstract_state_index()
        total_reward = 0
        
        for step in range(max_steps):
            valid_actions = game.get_valid_actions()
            action = self.choose_action(state_index, valid_actions)
            reward, next_state, done = game.step(action)
            next_state_index = game.get_current_abstract_state_index()
            
            # Add reward shaping
            shaped_reward = reward + self.add_reward_shaping(state_index, action, game)
            
            self.update_q_value(state_index, action, shaped_reward, next_state_index)
            
            total_reward += reward  # Track actual reward, not shaped
            state_index = next_state_index
            
            if done:
                break
        
        self.decay_epsilon()
        self.decay_learning_rate()
        self.episode_rewards.append(total_reward)
        self.recent_performance.append(total_reward)
        
        # Keep recent performance window manageable
        if len(self.recent_performance) > self.performance_window:
            self.recent_performance.pop(0)
        
        self.episode_count += 1
        
        return total_reward
    
    def train(self, game, num_episodes: int, max_steps_per_episode: int = 20, 
              verbose: bool = True, report_interval: int = 1000) -> List[float]:
        
        print(f"Starting enhanced Q-learning training for {num_episodes} episodes...")
        print(f"Enhanced features: Double Q-Learning={self.use_double_q}, Adaptive Epsilon={self.adaptive_epsilon}")
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
            if (episode + 1) % 100 == 0:
                recent_avg = np.mean(self.episode_rewards[-100:])
                recent_q_change = np.mean(self.q_value_changes[-100:]) if self.q_value_changes else 0
                pbar.set_postfix({
                    'avg_reward': f'{recent_avg:.2f}',
                    'epsilon': f'{self.epsilon:.3f}',
                    'lr': f'{self.learning_rate:.3f}',
                    'q_change': f'{recent_q_change:.4f}'
                })
            
            if verbose and (episode + 1) % report_interval == 0:
                avg_reward = np.mean(self.episode_rewards[-report_interval:])
                reward_std = np.std(self.episode_rewards[-report_interval:])
                best_reward = max(self.episode_rewards[-report_interval:])
                improvement_rate = self._calculate_improvement_rate()
                
                tqdm.write(f"\nEpisode {episode + 1}/{num_episodes}")
                tqdm.write(f"  Average reward (last {report_interval}): {avg_reward:.3f} ± {reward_std:.3f}")
                tqdm.write(f"  Best reward (last {report_interval}): {best_reward:.2f}")
                tqdm.write(f"  Current epsilon: {self.epsilon:.4f}")
                tqdm.write(f"  Current learning rate: {self.learning_rate:.4f}")
                tqdm.write(f"  Improvement rate: {improvement_rate:.4f}")
                
                if self.use_double_q:
                    sparsity_a = np.count_nonzero(self.q_table_a) / self.q_table_a.size * 100
                    sparsity_b = np.count_nonzero(self.q_table_b) / self.q_table_b.size * 100
                    tqdm.write(f"  Q-table sparsity: A={sparsity_a:.1f}%, B={sparsity_b:.1f}%")
                else:
                    sparsity = np.count_nonzero(self.q_table) / self.q_table.size * 100
                    tqdm.write(f"  Q-table sparsity: {sparsity:.1f}%")
        
        pbar.close()
        
        print(f"\nEnhanced training completed!")
        print(f"Final epsilon: {self.epsilon:.4f}")
        print(f"Final learning rate: {self.learning_rate:.4f}")
        print(f"Average reward (last 1000 episodes): {np.mean(self.episode_rewards[-1000:]):.3f}")
        print(f"Best reward achieved: {max(self.episode_rewards):.2f}")
        
        return self.episode_rewards
    
    def get_policy(self) -> np.ndarray:
        if self.use_double_q:
            # Use average of both Q-tables for final policy
            combined_q = (self.q_table_a + self.q_table_b) / 2
            return np.argmax(combined_q, axis=1)
        else:
            return np.argmax(self.q_table, axis=1)
    
    def evaluate_policy(self, game, num_episodes: int = 1000) -> Tuple[float, float]:
        old_epsilon = self.epsilon
        self.epsilon = 0.0  # No exploration during evaluation
        
        rewards = []
        
        # Create progress bar for evaluation
        pbar = tqdm(range(num_episodes), desc="Evaluating Enhanced Policy", leave=False)
        
        for _ in pbar:
            state = game.reset()
            total_reward = 0
            
            for _ in range(20):  # Max steps per episode
                state_index = game.get_current_abstract_state_index()
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
    
    def save_model(self, filename: str):
        if self.use_double_q:
            np.savez(filename, 
                     q_table_a=self.q_table_a,
                     q_table_b=self.q_table_b,
                     episode_rewards=np.array(self.episode_rewards),
                     epsilon=self.epsilon,
                     learning_rate=self.learning_rate,
                     episode_count=self.episode_count,
                     use_double_q=self.use_double_q)
        else:
            np.savez(filename, 
                     q_table=self.q_table,
                     episode_rewards=np.array(self.episode_rewards),
                     epsilon=self.epsilon,
                     learning_rate=self.learning_rate,
                     episode_count=self.episode_count,
                     use_double_q=self.use_double_q)
        print(f"Enhanced model saved to {filename}")
    
    def load_model(self, filename: str):
        data = np.load(filename)
        
        if data['use_double_q']:
            self.q_table_a = data['q_table_a']
            self.q_table_b = data['q_table_b']
        else:
            self.q_table = data['q_table']
            
        self.episode_rewards = data['episode_rewards'].tolist()
        self.epsilon = float(data['epsilon'])
        self.learning_rate = float(data['learning_rate'])
        self.episode_count = int(data['episode_count'])
        self.use_double_q = bool(data['use_double_q'])
        print(f"Enhanced model loaded from {filename}")