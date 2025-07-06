import numpy as np
import random
import heapq
from typing import Tuple, List, Dict, NamedTuple
from tqdm import tqdm
from collections import deque

class Experience(NamedTuple):
    """Experience tuple for replay buffer."""
    state: int
    action: int
    reward: float
    next_state: int
    done: bool
    priority: float

class PrioritizedReplayBuffer:
    """Prioritized Experience Replay buffer implementation."""
    
    def __init__(self, capacity: int, alpha: float = 0.6, beta: float = 0.4, beta_increment: float = 0.001):
        self.capacity = capacity
        self.alpha = alpha  # Priority exponent
        self.beta = beta    # Importance sampling exponent
        self.beta_increment = beta_increment
        self.buffer = []
        self.priorities = []
        self.position = 0
        
    def add(self, experience: Experience):
        """Add experience to buffer with priority."""
        priority = experience.priority if experience.priority > 0 else 1.0
        
        if len(self.buffer) < self.capacity:
            self.buffer.append(experience)
            self.priorities.append(priority)
        else:
            self.buffer[self.position] = experience
            self.priorities[self.position] = priority
            
        self.position = (self.position + 1) % self.capacity
    
    def sample(self, batch_size: int) -> Tuple[List[Experience], List[int], np.ndarray]:
        """Sample experiences based on priorities."""
        if len(self.buffer) < batch_size:
            batch_size = len(self.buffer)
        
        # Calculate sampling probabilities
        priorities_array = np.array(self.priorities[:len(self.buffer)])
        probs = priorities_array ** self.alpha
        probs = probs / probs.sum()
        
        # Sample indices
        indices = np.random.choice(len(self.buffer), batch_size, p=probs, replace=False)
        
        # Get experiences
        experiences = [self.buffer[i] for i in indices]
        
        # Calculate importance sampling weights
        weights = (len(self.buffer) * probs[indices]) ** (-self.beta)
        weights = weights / weights.max()  # Normalize
        
        # Update beta
        self.beta = min(1.0, self.beta + self.beta_increment)
        
        return experiences, indices, weights
    
    def update_priorities(self, indices: List[int], td_errors: List[float]):
        """Update priorities based on TD errors."""
        for idx, td_error in zip(indices, td_errors):
            if idx < len(self.priorities):
                self.priorities[idx] = abs(td_error) + 1e-6  # Small epsilon to avoid zero priority
    
    def __len__(self):
        return len(self.buffer)

class AdvancedQLearningAgent:
    """Advanced Q-Learning with UCB exploration, Prioritized Replay, and Multi-step learning."""
    
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
                 use_ucb: bool = True,
                 ucb_c: float = 2.0,
                 use_replay: bool = True,
                 replay_capacity: int = 50000,
                 batch_size: int = 32,
                 n_step: int = 3):
        
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
        self.use_ucb = use_ucb
        self.ucb_c = ucb_c
        self.use_replay = use_replay
        self.batch_size = batch_size
        self.n_step = n_step
        
        # Q-tables
        if use_double_q:
            self.q_table_a = np.zeros((state_space_size, action_space_size))
            self.q_table_b = np.zeros((state_space_size, action_space_size))
        else:
            self.q_table = np.zeros((state_space_size, action_space_size))
        
        # UCB exploration tracking
        if use_ucb:
            self.action_counts = np.zeros((state_space_size, action_space_size))
            self.state_visits = np.zeros(state_space_size)
        
        # Experience Replay
        if use_replay:
            self.replay_buffer = PrioritizedReplayBuffer(replay_capacity)
            self.n_step_buffer = deque(maxlen=n_step)
        
        # Tracking
        self.episode_rewards = []
        self.episode_count = 0
        self.recent_performance = []
        self.performance_window = 1000
        self.q_value_changes = []
        self.total_steps = 0
        
    def choose_action(self, state_index: int, valid_actions: List[int] = None) -> int:
        """Choose action using UCB or epsilon-greedy."""
        if valid_actions is None:
            valid_actions = list(range(self.action_space_size))
        
        # UCB exploration
        if self.use_ucb and self.total_steps > 100:  # Use UCB after initial exploration
            return self._choose_ucb_action(state_index, valid_actions)
        
        # Epsilon-greedy fallback
        if random.random() < self.epsilon:
            return random.choice(valid_actions)
        else:
            if self.use_double_q:
                combined_q = (self.q_table_a[state_index] + self.q_table_b[state_index]) / 2
            else:
                combined_q = self.q_table[state_index]
            
            masked_q_values = {action: combined_q[action] for action in valid_actions}
            return max(masked_q_values, key=masked_q_values.get)
    
    def _choose_ucb_action(self, state_index: int, valid_actions: List[int]) -> int:
        """Choose action using Upper Confidence Bound."""
        if self.use_double_q:
            combined_q = (self.q_table_a[state_index] + self.q_table_b[state_index]) / 2
        else:
            combined_q = self.q_table[state_index]
        
        ucb_values = {}
        total_visits = max(1, self.state_visits[state_index])
        
        for action in valid_actions:
            action_count = max(1, self.action_counts[state_index, action])
            
            # UCB formula: Q(s,a) + c * sqrt(ln(N(s)) / N(s,a))
            confidence_bonus = self.ucb_c * np.sqrt(np.log(total_visits) / action_count)
            ucb_value = combined_q[action] + confidence_bonus
            ucb_values[action] = ucb_value
        
        return max(ucb_values, key=ucb_values.get)
    
    def update_q_value(self, state_index: int, action: int, reward: float, next_state_index: int, done: bool = False):
        """Update Q-values with optional experience replay."""
        
        # Store experience for replay
        if self.use_replay:
            # Calculate initial TD error for priority
            if self.use_double_q:
                current_q = (self.q_table_a[state_index, action] + self.q_table_b[state_index, action]) / 2
                next_q = np.max((self.q_table_a[next_state_index] + self.q_table_b[next_state_index]) / 2)
            else:
                current_q = self.q_table[state_index, action]
                next_q = np.max(self.q_table[next_state_index])
            
            target = reward + self.discount_factor * next_q * (not done)
            td_error = abs(target - current_q)
            
            experience = Experience(state_index, action, reward, next_state_index, done, td_error)
            
            # Add to n-step buffer
            self.n_step_buffer.append(experience)
            
            # If we have enough steps, compute n-step return and add to replay buffer
            if len(self.n_step_buffer) == self.n_step:
                n_step_exp = self._compute_n_step_experience()
                self.replay_buffer.add(n_step_exp)
        
        # Direct Q-learning update
        self._update_q_direct(state_index, action, reward, next_state_index, done)
        
        # Experience replay training
        if self.use_replay and len(self.replay_buffer) >= self.batch_size:
            self._replay_training()
        
        # Update UCB counts
        if self.use_ucb:
            self.action_counts[state_index, action] += 1
            self.state_visits[state_index] += 1
    
    def _compute_n_step_experience(self) -> Experience:
        """Compute n-step return from buffer."""
        if len(self.n_step_buffer) == 0:
            return None
        
        # Get first experience
        first_exp = self.n_step_buffer[0]
        
        # Compute n-step return
        n_step_return = 0
        gamma = 1
        done = False
        
        for exp in self.n_step_buffer:
            n_step_return += gamma * exp.reward
            gamma *= self.discount_factor
            if exp.done:
                done = True
                break
        
        # Final state is the last state in buffer
        final_state = self.n_step_buffer[-1].next_state
        
        # Calculate priority (TD error)
        if self.use_double_q:
            current_q = (self.q_table_a[first_exp.state, first_exp.action] + 
                        self.q_table_b[first_exp.state, first_exp.action]) / 2
            if not done:
                next_q = np.max((self.q_table_a[final_state] + self.q_table_b[final_state]) / 2)
                target = n_step_return + gamma * next_q
            else:
                target = n_step_return
        else:
            current_q = self.q_table[first_exp.state, first_exp.action]
            if not done:
                next_q = np.max(self.q_table[final_state])
                target = n_step_return + gamma * next_q
            else:
                target = n_step_return
        
        td_error = abs(target - current_q)
        
        return Experience(first_exp.state, first_exp.action, n_step_return, 
                         final_state, done, td_error)
    
    def _update_q_direct(self, state_index: int, action: int, reward: float, next_state_index: int, done: bool):
        """Direct Q-learning update."""
        if self.use_double_q:
            # Double Q-Learning update
            if random.random() < 0.5:
                old_q = self.q_table_a[state_index, action]
                best_next_action = np.argmax(self.q_table_a[next_state_index])
                max_next_q = self.q_table_b[next_state_index, best_next_action] * (not done)
                
                new_q = old_q + self.learning_rate * (
                    reward + self.discount_factor * max_next_q - old_q
                )
                self.q_table_a[state_index, action] = new_q
            else:
                old_q = self.q_table_b[state_index, action]
                best_next_action = np.argmax(self.q_table_b[next_state_index])
                max_next_q = self.q_table_a[next_state_index, best_next_action] * (not done)
                
                new_q = old_q + self.learning_rate * (
                    reward + self.discount_factor * max_next_q - old_q
                )
                self.q_table_b[state_index, action] = new_q
        else:
            old_q = self.q_table[state_index, action]
            max_next_q = np.max(self.q_table[next_state_index]) * (not done)
            
            new_q = old_q + self.learning_rate * (
                reward + self.discount_factor * max_next_q - old_q
            )
            self.q_table[state_index, action] = new_q
        
        # Track Q-value changes
        q_change = abs(new_q - old_q)
        self.q_value_changes.append(q_change)
        if len(self.q_value_changes) > 1000:
            self.q_value_changes.pop(0)
    
    def _replay_training(self):
        """Train on batch of experiences from replay buffer."""
        if len(self.replay_buffer) < self.batch_size:
            return
        
        # Sample experiences
        experiences, indices, weights = self.replay_buffer.sample(self.batch_size)
        
        td_errors = []
        
        for i, exp in enumerate(experiences):
            # Compute target value
            if self.use_double_q:
                if random.random() < 0.5:
                    current_q = self.q_table_a[exp.state, exp.action]
                    if not exp.done:
                        best_action = np.argmax(self.q_table_a[exp.next_state])
                        target_q = self.q_table_b[exp.next_state, best_action]
                    else:
                        target_q = 0
                    
                    target = exp.reward + self.discount_factor * target_q
                    td_error = target - current_q
                    
                    # Weighted update
                    new_q = current_q + self.learning_rate * weights[i] * td_error
                    self.q_table_a[exp.state, exp.action] = new_q
                else:
                    current_q = self.q_table_b[exp.state, exp.action]
                    if not exp.done:
                        best_action = np.argmax(self.q_table_b[exp.next_state])
                        target_q = self.q_table_a[exp.next_state, best_action]
                    else:
                        target_q = 0
                    
                    target = exp.reward + self.discount_factor * target_q
                    td_error = target - current_q
                    
                    # Weighted update
                    new_q = current_q + self.learning_rate * weights[i] * td_error
                    self.q_table_b[exp.state, exp.action] = new_q
            else:
                current_q = self.q_table[exp.state, exp.action]
                if not exp.done:
                    target_q = np.max(self.q_table[exp.next_state])
                else:
                    target_q = 0
                
                target = exp.reward + self.discount_factor * target_q
                td_error = target - current_q
                
                # Weighted update
                new_q = current_q + self.learning_rate * weights[i] * td_error
                self.q_table[exp.state, exp.action] = new_q
            
            td_errors.append(abs(td_error))
        
        # Update priorities
        self.replay_buffer.update_priorities(indices, td_errors)
    
    def decay_epsilon(self):
        """Adaptive epsilon decay."""
        if len(self.recent_performance) > 100:
            recent_improvement = self._calculate_improvement_rate()
            
            if recent_improvement > 0.1:
                decay_rate = self.epsilon_decay * 0.99
            elif recent_improvement < -0.05:
                decay_rate = self.epsilon_decay * 1.001
            else:
                decay_rate = self.epsilon_decay
            
            if self.epsilon > self.epsilon_min:
                self.epsilon *= decay_rate
        else:
            if self.epsilon > self.epsilon_min:
                self.epsilon *= self.epsilon_decay
    
    def _calculate_improvement_rate(self) -> float:
        """Calculate recent performance improvement rate."""
        if len(self.recent_performance) < 200:
            return 0.0
        
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
        
        counts = game_state.get_current_counts()
        snake_count = counts[3]
        net_count = counts[4]
        clover_count = counts[6]
        
        if action == 8:  # Cash out
            if clover_count > 0:
                bonus += 0.1 * clover_count
            if snake_count > 0 and net_count == 0:
                bonus += 0.2
        elif action == 3:  # Respin snake
            if net_count == 0:
                bonus += 0.1
        elif action == 0:  # Respin blank
            if clover_count > 0 or (snake_count > 0 and net_count > 0):
                bonus += 0.05
        
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
            
            self.update_q_value(state_index, action, shaped_reward, next_state_index, done)
            
            total_reward += reward
            state_index = next_state_index
            self.total_steps += 1
            
            if done:
                break
        
        # Handle end of episode for n-step learning
        if self.use_replay:
            # Add remaining experiences from n-step buffer
            while len(self.n_step_buffer) > 0:
                exp = self.n_step_buffer.popleft()
                self.replay_buffer.add(exp)
        
        self.decay_epsilon()
        self.decay_learning_rate()
        self.episode_rewards.append(total_reward)
        self.recent_performance.append(total_reward)
        
        if len(self.recent_performance) > self.performance_window:
            self.recent_performance.pop(0)
        
        self.episode_count += 1
        
        return total_reward
    
    def train(self, game, num_episodes: int, max_steps_per_episode: int = 20, 
              verbose: bool = True, report_interval: int = 1000) -> List[float]:
        
        print(f"Starting advanced Q-learning training for {num_episodes} episodes...")
        print(f"Advanced features:")
        print(f"  Double Q-Learning: {self.use_double_q}")
        print(f"  UCB Exploration: {self.use_ucb}")
        print(f"  Prioritized Replay: {self.use_replay}")
        print(f"  N-step Learning: {self.n_step}")
        print(f"Initial epsilon: {self.epsilon:.3f}")
        print(f"Learning rate: {self.learning_rate}")
        print("-" * 50)
        
        pbar = tqdm(range(num_episodes), desc="Advanced Training", 
                   bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')
        
        for episode in pbar:
            episode_reward = self.train_episode(game, max_steps_per_episode)
            
            if (episode + 1) % 100 == 0:
                recent_avg = np.mean(self.episode_rewards[-100:])
                recent_q_change = np.mean(self.q_value_changes[-100:]) if self.q_value_changes else 0
                
                postfix = {
                    'avg_reward': f'{recent_avg:.2f}',
                    'epsilon': f'{self.epsilon:.3f}',
                    'lr': f'{self.learning_rate:.3f}',
                    'q_change': f'{recent_q_change:.4f}'
                }
                
                if self.use_replay:
                    postfix['replay_size'] = len(self.replay_buffer)
                
                pbar.set_postfix(postfix)
            
            if verbose and (episode + 1) % report_interval == 0:
                avg_reward = np.mean(self.episode_rewards[-report_interval:])
                reward_std = np.std(self.episode_rewards[-report_interval:])
                best_reward = max(self.episode_rewards[-report_interval:])
                improvement_rate = self._calculate_improvement_rate()
                
                tqdm.write(f"\nEpisode {episode + 1}/{num_episodes}")
                tqdm.write(f"  Average reward (last {report_interval}): {avg_reward:.3f} ± {reward_std:.3f}")
                tqdm.write(f"  Best reward: {best_reward:.2f}")
                tqdm.write(f"  Current epsilon: {self.epsilon:.4f}")
                tqdm.write(f"  Learning rate: {self.learning_rate:.4f}")
                tqdm.write(f"  Improvement rate: {improvement_rate:.4f}")
                
                if self.use_replay:
                    tqdm.write(f"  Replay buffer size: {len(self.replay_buffer)}")
                
                if self.use_double_q:
                    sparsity_a = np.count_nonzero(self.q_table_a) / self.q_table_a.size * 100
                    sparsity_b = np.count_nonzero(self.q_table_b) / self.q_table_b.size * 100
                    tqdm.write(f"  Q-table sparsity: A={sparsity_a:.1f}%, B={sparsity_b:.1f}%")
        
        pbar.close()
        
        print(f"\nAdvanced training completed!")
        print(f"Final epsilon: {self.epsilon:.4f}")
        print(f"Final learning rate: {self.learning_rate:.4f}")
        print(f"Average reward (last 1000): {np.mean(self.episode_rewards[-1000:]):.3f}")
        print(f"Best reward achieved: {max(self.episode_rewards):.2f}")
        
        if self.use_replay:
            print(f"Final replay buffer size: {len(self.replay_buffer)}")
        
        return self.episode_rewards
    
    def get_policy(self) -> np.ndarray:
        if self.use_double_q:
            combined_q = (self.q_table_a + self.q_table_b) / 2
            return np.argmax(combined_q, axis=1)
        else:
            return np.argmax(self.q_table, axis=1)
    
    def evaluate_policy(self, game, num_episodes: int = 1000) -> Tuple[float, float]:
        old_epsilon = self.epsilon
        self.epsilon = 0.0
        
        rewards = []
        pbar = tqdm(range(num_episodes), desc="Evaluating Advanced Policy", leave=False)
        
        for _ in pbar:
            state = game.reset()
            total_reward = 0
            
            for _ in range(20):
                state_index = game.get_current_abstract_state_index()
                valid_actions = game.get_valid_actions()
                action = self.choose_action(state_index, valid_actions)
                reward, next_state, done = game.step(action)
                total_reward += reward
                
                if done:
                    break
            
            rewards.append(total_reward)
            
            if len(rewards) % 100 == 0:
                current_avg = np.mean(rewards)
                pbar.set_postfix({'avg_reward': f'{current_avg:.3f}'})
        
        pbar.close()
        self.epsilon = old_epsilon
        return np.mean(rewards), np.std(rewards)
    
    def save_model(self, filename: str):
        save_dict = {
            'episode_rewards': np.array(self.episode_rewards),
            'epsilon': self.epsilon,
            'learning_rate': self.learning_rate,
            'episode_count': self.episode_count,
            'use_double_q': self.use_double_q,
            'use_ucb': self.use_ucb,
            'use_replay': self.use_replay,
            'n_step': self.n_step,
            'total_steps': self.total_steps
        }
        
        if self.use_double_q:
            save_dict['q_table_a'] = self.q_table_a
            save_dict['q_table_b'] = self.q_table_b
        else:
            save_dict['q_table'] = self.q_table
            
        if self.use_ucb:
            save_dict['action_counts'] = self.action_counts
            save_dict['state_visits'] = self.state_visits
        
        np.savez(filename, **save_dict)
        print(f"Advanced model saved to {filename}")
    
    def load_model(self, filename: str):
        data = np.load(filename)
        
        self.episode_rewards = data['episode_rewards'].tolist()
        self.epsilon = float(data['epsilon'])
        self.learning_rate = float(data['learning_rate'])
        self.episode_count = int(data['episode_count'])
        self.use_double_q = bool(data['use_double_q'])
        self.use_ucb = bool(data['use_ucb'])
        self.use_replay = bool(data['use_replay'])
        self.n_step = int(data['n_step'])
        self.total_steps = int(data['total_steps'])
        
        if self.use_double_q:
            self.q_table_a = data['q_table_a']
            self.q_table_b = data['q_table_b']
        else:
            self.q_table = data['q_table']
            
        if self.use_ucb:
            self.action_counts = data['action_counts']
            self.state_visits = data['state_visits']
        
        print(f"Advanced model loaded from {filename}")