import numpy as np
from typing import List, Tuple, Dict
from advanced_q_learning import AdvancedQLearningAgent
from enhanced_q_learning import EnhancedQLearningAgent
from q_learning import QLearningAgent
from tqdm import tqdm

class EnsembleQLearningAgent:
    """Ensemble of Q-learning agents with different configurations."""
    
    def __init__(self, state_space_size: int, action_space_size: int):
        self.state_space_size = state_space_size
        self.action_space_size = action_space_size
        self.agents = []
        self.agent_weights = []
        self.episode_rewards = []
        self.episode_count = 0
        
        self._create_ensemble()
    
    def _create_ensemble(self):
        """Create diverse ensemble of agents."""
        
        # Agent 1: Conservative baseline
        agent1 = QLearningAgent(
            state_space_size=self.state_space_size,
            action_space_size=self.action_space_size,
            learning_rate=0.1,
            discount_factor=0.95,
            epsilon=1.0,
            epsilon_decay=0.999,
            epsilon_min=0.05
        )
        
        # Agent 2: Aggressive learner
        agent2 = EnhancedQLearningAgent(
            state_space_size=self.state_space_size,
            action_space_size=self.action_space_size,
            learning_rate=0.4,
            discount_factor=0.99,
            epsilon=1.0,
            epsilon_decay=0.9995,
            epsilon_min=0.01,
            use_double_q=True,
            adaptive_epsilon=True
        )
        
        # Agent 3: UCB explorer
        agent3 = AdvancedQLearningAgent(
            state_space_size=self.state_space_size,
            action_space_size=self.action_space_size,
            learning_rate=0.25,
            discount_factor=0.99,
            epsilon=0.5,
            epsilon_decay=0.995,
            epsilon_min=0.1,
            use_double_q=True,
            use_ucb=True,
            ucb_c=1.5,
            use_replay=False,  # Focus on UCB
            n_step=1
        )
        
        # Agent 4: Experience replay specialist
        agent4 = AdvancedQLearningAgent(
            state_space_size=self.state_space_size,
            action_space_size=self.action_space_size,
            learning_rate=0.2,
            discount_factor=0.99,
            epsilon=1.0,
            epsilon_decay=0.9995,
            epsilon_min=0.02,
            use_double_q=True,
            use_ucb=False,
            use_replay=True,
            replay_capacity=10000,
            batch_size=64,
            n_step=5
        )
        
        # Agent 5: Multi-step specialist
        agent5 = AdvancedQLearningAgent(
            state_space_size=self.state_space_size,
            action_space_size=self.action_space_size,
            learning_rate=0.3,
            discount_factor=0.98,
            epsilon=1.0,
            epsilon_decay=0.9990,
            epsilon_min=0.01,
            use_double_q=True,
            use_ucb=False,
            use_replay=True,
            replay_capacity=5000,
            batch_size=32,
            n_step=7
        )
        
        self.agents = [agent1, agent2, agent3, agent4, agent5]
        self.agent_weights = [1.0] * len(self.agents)  # Equal weights initially
        
        # Agent descriptions for reporting
        self.agent_descriptions = [
            "Conservative Baseline",
            "Aggressive Enhanced",
            "UCB Explorer", 
            "Replay Specialist",
            "Multi-step Specialist"
        ]
    
    def choose_action(self, state_index: int, valid_actions: List[int] = None) -> int:
        """Choose action using weighted voting from ensemble."""
        if valid_actions is None:
            valid_actions = list(range(self.action_space_size))
        
        # Get action preferences from each agent
        action_votes = {}
        for action in valid_actions:
            action_votes[action] = 0.0
        
        total_weight = sum(self.agent_weights)
        
        for i, agent in enumerate(self.agents):
            # Get agent's preferred action
            agent_action = agent.choose_action(state_index, valid_actions)
            
            # Weight the vote
            weight = self.agent_weights[i] / total_weight
            action_votes[agent_action] += weight
        
        # Choose action with highest weighted vote
        return max(action_votes, key=action_votes.get)
    
    def update_ensemble_weights(self, episode_rewards: List[float]):
        """Update agent weights based on recent performance."""
        if len(episode_rewards) < 100:
            return
        
        # Calculate recent performance for each agent
        recent_window = 100
        agent_performances = []
        
        for i, agent in enumerate(self.agents):
            if len(agent.episode_rewards) >= recent_window:
                recent_perf = np.mean(agent.episode_rewards[-recent_window:])
                agent_performances.append(recent_perf)
            else:
                agent_performances.append(np.mean(agent.episode_rewards) if agent.episode_rewards else 0)
        
        # Convert to positive weights (shift if needed)
        min_perf = min(agent_performances)
        if min_perf < 0:
            agent_performances = [p - min_perf + 1 for p in agent_performances]
        
        # Normalize weights with temperature for softmax-like behavior
        temperature = 2.0
        exp_perfs = [np.exp(p / temperature) for p in agent_performances]
        total_exp = sum(exp_perfs)
        
        if total_exp > 0:
            self.agent_weights = [exp_p / total_exp * len(self.agents) for exp_p in exp_perfs]
        else:
            self.agent_weights = [1.0] * len(self.agents)
    
    def train_episode(self, game, max_steps: int = 20) -> float:
        """Train all agents on the same episode."""
        # Reset game
        initial_state = game.reset()
        
        # Store episode trajectory for all agents
        episode_trajectory = []
        
        # Play episode
        state = initial_state
        state_index = game.get_current_abstract_state_index()
        total_reward = 0
        
        for step in range(max_steps):
            valid_actions = game.get_valid_actions()
            
            # Choose action using ensemble
            action = self.choose_action(state_index, valid_actions)
            
            # Take action
            reward, next_state, done = game.step(action)
            next_state_index = game.get_current_abstract_state_index()
            
            # Store step
            episode_trajectory.append({
                'state_index': state_index,
                'action': action,
                'reward': reward,
                'next_state_index': next_state_index,
                'done': done,
                'valid_actions': valid_actions.copy()
            })
            
            total_reward += reward
            state_index = next_state_index
            
            if done:
                break
        
        # Train all agents on the same trajectory
        for agent in self.agents:
            agent_reward = self._train_agent_on_trajectory(agent, game, episode_trajectory)
        
        # Update ensemble weights periodically
        if self.episode_count % 100 == 0:
            all_rewards = [agent.episode_rewards for agent in self.agents if agent.episode_rewards]
            if all_rewards:
                self.update_ensemble_weights(all_rewards[0])  # Use first agent's rewards as reference
        
        self.episode_rewards.append(total_reward)
        self.episode_count += 1
        
        return total_reward
    
    def _train_agent_on_trajectory(self, agent, game, trajectory: List[Dict]) -> float:
        """Train individual agent on episode trajectory."""
        total_reward = 0
        
        for i, step in enumerate(trajectory):
            # Reset game to replay the step
            game.current_state = game.index_to_state(step['state_index'])
            game.respins_used = step['state_index'] % (game.max_respins + 1)  # Simplified
            
            # Let agent choose its own action for learning
            agent_action = agent.choose_action(step['state_index'], step['valid_actions'])
            
            # Update agent with actual trajectory action and reward for consistency
            if hasattr(agent, 'add_reward_shaping'):
                shaped_reward = step['reward'] + agent.add_reward_shaping(step['state_index'], step['action'], game)
            else:
                shaped_reward = step['reward']
            
            # Update agent's Q-values
            if hasattr(agent, 'update_q_value'):
                agent.update_q_value(step['state_index'], step['action'], shaped_reward, 
                                   step['next_state_index'], step['done'])
            
            total_reward += step['reward']
        
        # Update agent's episode tracking
        if hasattr(agent, 'decay_epsilon'):
            agent.decay_epsilon()
        if hasattr(agent, 'decay_learning_rate'):
            agent.decay_learning_rate()
        
        agent.episode_rewards.append(total_reward)
        agent.episode_count += 1
        
        if hasattr(agent, 'recent_performance'):
            agent.recent_performance.append(total_reward)
            if len(agent.recent_performance) > getattr(agent, 'performance_window', 1000):
                agent.recent_performance.pop(0)
        
        return total_reward
    
    def train(self, game, num_episodes: int, max_steps_per_episode: int = 20, 
              verbose: bool = True, report_interval: int = 1000) -> List[float]:
        
        print(f"Starting ensemble Q-learning training for {num_episodes} episodes...")
        print(f"Ensemble composition:")
        for i, desc in enumerate(self.agent_descriptions):
            print(f"  Agent {i+1}: {desc}")
        print("-" * 50)
        
        pbar = tqdm(range(num_episodes), desc="Ensemble Training", 
                   bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]')
        
        for episode in pbar:
            episode_reward = self.train_episode(game, max_steps_per_episode)
            
            if (episode + 1) % 100 == 0:
                recent_avg = np.mean(self.episode_rewards[-100:])
                
                # Show current weights
                weight_str = ", ".join([f"{w:.2f}" for w in self.agent_weights])
                
                pbar.set_postfix({
                    'avg_reward': f'{recent_avg:.2f}',
                    'weights': f'[{weight_str}]'
                })
            
            if verbose and (episode + 1) % report_interval == 0:
                avg_reward = np.mean(self.episode_rewards[-report_interval:])
                reward_std = np.std(self.episode_rewards[-report_interval:])
                best_reward = max(self.episode_rewards[-report_interval:])
                
                tqdm.write(f"\nEpisode {episode + 1}/{num_episodes}")
                tqdm.write(f"  Ensemble avg reward (last {report_interval}): {avg_reward:.3f} ± {reward_std:.3f}")
                tqdm.write(f"  Best reward: {best_reward:.2f}")
                tqdm.write(f"  Current agent weights:")
                
                for i, (desc, weight) in enumerate(zip(self.agent_descriptions, self.agent_weights)):
                    agent_recent = np.mean(self.agents[i].episode_rewards[-100:]) if len(self.agents[i].episode_rewards) >= 100 else 0
                    tqdm.write(f"    {desc}: {weight:.3f} (recent: {agent_recent:.2f})")
        
        pbar.close()
        
        print(f"\nEnsemble training completed!")
        print(f"Final ensemble weights:")
        for i, (desc, weight) in enumerate(zip(self.agent_descriptions, self.agent_weights)):
            print(f"  {desc}: {weight:.3f}")
        print(f"Average reward (last 1000): {np.mean(self.episode_rewards[-1000:]):.3f}")
        print(f"Best reward achieved: {max(self.episode_rewards):.2f}")
        
        return self.episode_rewards
    
    def get_policy(self) -> np.ndarray:
        """Get ensemble policy by weighted voting."""
        ensemble_q = np.zeros((self.state_space_size, self.action_space_size))
        total_weight = sum(self.agent_weights)
        
        for i, agent in enumerate(self.agents):
            weight = self.agent_weights[i] / total_weight
            
            if hasattr(agent, 'q_table_a') and agent.use_double_q:
                agent_q = (agent.q_table_a + agent.q_table_b) / 2
            elif hasattr(agent, 'q_table'):
                agent_q = agent.q_table
            else:
                continue
            
            ensemble_q += weight * agent_q
        
        return np.argmax(ensemble_q, axis=1)
    
    def evaluate_policy(self, game, num_episodes: int = 1000) -> Tuple[float, float]:
        """Evaluate ensemble policy."""
        # Temporarily disable exploration for all agents
        old_epsilons = []
        for agent in self.agents:
            old_epsilons.append(agent.epsilon)
            agent.epsilon = 0.0
        
        rewards = []
        pbar = tqdm(range(num_episodes), desc="Evaluating Ensemble Policy", leave=False)
        
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
        
        # Restore original epsilons
        for agent, old_eps in zip(self.agents, old_epsilons):
            agent.epsilon = old_eps
        
        return np.mean(rewards), np.std(rewards)
    
    def save_model(self, filename: str):
        """Save ensemble model."""
        # Save individual agents
        for i, agent in enumerate(self.agents):
            agent_filename = f"{filename}_agent_{i}.npz"
            agent.save_model(agent_filename)
        
        # Save ensemble data
        ensemble_data = {
            'agent_weights': np.array(self.agent_weights),
            'episode_rewards': np.array(self.episode_rewards),
            'episode_count': self.episode_count,
            'agent_descriptions': self.agent_descriptions
        }
        
        np.savez(f"{filename}_ensemble.npz", **ensemble_data)
        print(f"Ensemble model saved with base filename: {filename}")
    
    def load_model(self, filename: str):
        """Load ensemble model."""
        # Load ensemble data
        ensemble_data = np.load(f"{filename}_ensemble.npz", allow_pickle=True)
        self.agent_weights = ensemble_data['agent_weights'].tolist()
        self.episode_rewards = ensemble_data['episode_rewards'].tolist()
        self.episode_count = int(ensemble_data['episode_count'])
        self.agent_descriptions = ensemble_data['agent_descriptions'].tolist()
        
        # Load individual agents
        for i, agent in enumerate(self.agents):
            agent_filename = f"{filename}_agent_{i}.npz"
            try:
                agent.load_model(agent_filename)
            except FileNotFoundError:
                print(f"Warning: Could not load agent {i} from {agent_filename}")
        
        print(f"Ensemble model loaded from base filename: {filename}")