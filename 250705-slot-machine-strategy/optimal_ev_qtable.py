import numpy as np
from typing import Tuple, List, Dict
from itertools import product
from tqdm import tqdm
from game_rules import SlotMachine

class OptimalEVQTable:
    """Pre-computed Q-table with optimal expected values for all states."""
    
    def __init__(self, game: SlotMachine):
        self.game = game
        self.symbols = game.symbols
        self.symbol_probabilities = game.symbol_probabilities
        self.num_symbols = len(self.symbols)
        self.max_respins = game.max_respins
        
        # Q-table: [state_index][action] = expected_value
        self.q_table = None
        self.state_space_size = self._calculate_state_space_size()
        self.action_space_size = game.get_action_space_size()
        
        # Memoization for recursive calculations
        self.memo = {}
        
    def _calculate_state_space_size(self) -> int:
        """Calculate the size of the abstract state space."""
        # We use symbol counts (0-4 for each of 8 symbols) + respins_used (0-5)
        # This gives us all possible combinations of symbol counts that sum to 4
        return 500 * (self.max_respins + 1)  # Using same abstraction as game
    
    def _counts_to_state(self, symbol_counts: List[int], respins_used: int) -> Tuple[str, ...]:
        """Convert symbol counts to a representative state tuple."""
        state = []
        for symbol_idx, count in enumerate(symbol_counts):
            state.extend([self.symbols[symbol_idx]] * count)
        
        # Pad with blanks if needed
        while len(state) < 4:
            state.append('Blank')
        
        return tuple(state[:4])
    
    def _generate_all_possible_states(self) -> List[Tuple[List[int], int]]:
        """Generate all possible (symbol_counts, respins_used) combinations."""
        states = []
        
        # Generate all ways to distribute 4 symbols among 8 types
        for counts in product(range(5), repeat=8):  # 0-4 for each symbol type
            if sum(counts) == 4:  # Must sum to 4 (4 reels)
                for respins in range(self.max_respins + 1):
                    states.append((list(counts), respins))
        
        return states
    
    def _calculate_payout_from_counts(self, symbol_counts: List[int]) -> int:
        """Calculate payout from symbol counts."""
        state = self._counts_to_state(symbol_counts, 0)
        return self.game._calculate_payout(state, 1)
    
    def _calculate_expected_value_recursive(self, symbol_counts: List[int], respins_used: int) -> float:
        """Calculate expected value for a state using dynamic programming."""
        # Create a hashable key for memoization
        state_key = (tuple(symbol_counts), respins_used)
        
        if state_key in self.memo:
            return self.memo[state_key]
        
        # Base case: no respins left, must cash out
        if respins_used >= self.max_respins:
            cash_out_value = self._calculate_payout_from_counts(symbol_counts)
            self.memo[state_key] = cash_out_value
            return cash_out_value
        
        # Calculate cash-out value
        cash_out_value = self._calculate_payout_from_counts(symbol_counts)
        
        # Calculate expected value for each respin action
        best_value = cash_out_value
        
        # Try respinning each symbol type that's present
        for symbol_idx in range(self.num_symbols):
            if symbol_counts[symbol_idx] > 0:  # This symbol is present
                # Calculate expected value of respinning this symbol type
                respin_value = self._calculate_respin_expected_value(symbol_counts, symbol_idx, respins_used)
                best_value = max(best_value, respin_value)
        
        self.memo[state_key] = best_value
        return best_value
    
    def _calculate_respin_expected_value(self, symbol_counts: List[int], symbol_to_respin: int, respins_used: int) -> float:
        """Calculate expected value of respinning a specific symbol type."""
        expected_value = 0.0
        
        # For each possible new symbol that could replace the respun symbol
        for new_symbol_idx, prob in enumerate(self.symbol_probabilities):
            new_counts = symbol_counts.copy()
            new_counts[symbol_to_respin] -= 1  # Remove one of the old symbol
            new_counts[new_symbol_idx] += 1    # Add one of the new symbol
            
            # Calculate expected value of the new state
            new_state_value = self._calculate_expected_value_recursive(new_counts, respins_used + 1)
            expected_value += prob * new_state_value
        
        # Subtract the cost of respinning (1 Gold)
        return expected_value - 1
    
    def _get_best_action_for_state(self, symbol_counts: List[int], respins_used: int) -> Tuple[int, float]:
        """Get the best action and its expected value for a given state."""
        if respins_used >= self.max_respins:
            return 8, self._calculate_payout_from_counts(symbol_counts)  # Must cash out
        
        # Calculate cash-out value
        cash_out_value = self._calculate_payout_from_counts(symbol_counts)
        best_action = 8
        best_value = cash_out_value
        
        # Check each respin action
        for symbol_idx in range(self.num_symbols):
            if symbol_counts[symbol_idx] > 0:  # This symbol is present
                respin_value = self._calculate_respin_expected_value(symbol_counts, symbol_idx, respins_used)
                if respin_value > best_value:
                    best_value = respin_value
                    best_action = symbol_idx  # Action index matches symbol index
        
        return best_action, best_value
    
    def build_optimal_qtable(self) -> np.ndarray:
        """Build the optimal Q-table by calculating expected values for all states."""
        print("Building Optimal Expected Value Q-Table...")
        
        # Generate all possible states
        all_states = self._generate_all_possible_states()
        print(f"Total states to compute: {len(all_states)}")
        
        # Initialize Q-table
        q_table = np.zeros((self.state_space_size, self.action_space_size))
        
        # Calculate optimal values for each state
        for i, (symbol_counts, respins_used) in enumerate(tqdm(all_states, desc="Computing EV")):
            # Get state index using game's abstraction method
            state_tuple = self._counts_to_state(symbol_counts, respins_used)
            
            # Set game state temporarily to get proper index
            old_state = self.game.current_state
            old_respins = self.game.respins_used
            self.game.current_state = state_tuple
            self.game.respins_used = respins_used
            
            try:
                state_index = self.game.get_current_abstract_state_index()
                
                # Calculate Q-values for each action
                for action in range(self.action_space_size):
                    if action == 8:  # Cash out
                        q_table[state_index, action] = self._calculate_payout_from_counts(symbol_counts)
                    else:  # Respin actions
                        if symbol_counts[action] > 0:  # Symbol is present
                            q_table[state_index, action] = self._calculate_respin_expected_value(
                                symbol_counts, action, respins_used
                            )
                        else:  # Symbol not present - invalid action
                            q_table[state_index, action] = float('-inf')
                
            except Exception as e:
                print(f"Error processing state {i}: {e}")
                continue
            finally:
                # Restore game state
                self.game.current_state = old_state
                self.game.respins_used = old_respins
        
        self.q_table = q_table
        print("Q-table construction completed!")
        return q_table
    
    def get_optimal_action(self, symbol_counts: List[int], respins_used: int) -> Tuple[int, float]:
        """Get the optimal action for a given state using the pre-computed Q-table."""
        if self.q_table is None:
            raise ValueError("Q-table not built yet. Call build_optimal_qtable() first.")
        
        # Get state index
        state_tuple = self._counts_to_state(symbol_counts, respins_used)
        old_state = self.game.current_state
        old_respins = self.game.respins_used
        self.game.current_state = state_tuple
        self.game.respins_used = respins_used
        
        try:
            state_index = self.game.get_current_abstract_state_index()
            
            # Get valid actions (only those corresponding to symbols present)
            valid_actions = []
            for action in range(self.action_space_size):
                if action == 8:  # Cash out is always valid
                    valid_actions.append(action)
                elif action < len(symbol_counts) and symbol_counts[action] > 0:
                    valid_actions.append(action)
            
            # Find best action among valid ones
            best_action = 8
            best_value = self.q_table[state_index, 8]  # Cash out value
            
            for action in valid_actions:
                if self.q_table[state_index, action] > best_value:
                    best_value = self.q_table[state_index, action]
                    best_action = action
            
            return best_action, best_value
            
        finally:
            self.game.current_state = old_state
            self.game.respins_used = old_respins
    
    def analyze_qtable(self) -> Dict[str, any]:
        """Analyze the completed Q-table for insights."""
        if self.q_table is None:
            return {'error': 'Q-table not built yet'}
        
        # Count positive EV actions
        positive_ev_actions = np.sum(self.q_table > 0)
        total_actions = np.sum(self.q_table > float('-inf'))
        
        # Find states with positive EV respin actions
        cash_out_values = self.q_table[:, 8]  # Cash out column
        respin_values = self.q_table[:, :8]   # Respin actions
        
        # For each state, check if any respin action is better than cash out
        states_with_positive_ev = 0
        for state_idx in range(self.q_table.shape[0]):
            cash_out_val = cash_out_values[state_idx]
            best_respin = np.max(respin_values[state_idx])
            if best_respin > cash_out_val:
                states_with_positive_ev += 1
        
        # Calculate value ranges
        valid_values = self.q_table[self.q_table > float('-inf')]
        
        analysis = {
            'total_states': self.q_table.shape[0],
            'total_actions': self.q_table.shape[1],
            'positive_ev_actions': positive_ev_actions,
            'total_valid_actions': total_actions,
            'positive_ev_percentage': (positive_ev_actions / total_actions) * 100,
            'states_with_positive_ev_respins': states_with_positive_ev,
            'percentage_states_with_positive_ev': (states_with_positive_ev / self.q_table.shape[0]) * 100,
            'value_statistics': {
                'min_value': np.min(valid_values),
                'max_value': np.max(valid_values),
                'mean_value': np.mean(valid_values),
                'std_value': np.std(valid_values)
            }
        }
        
        return analysis
    
    def save_qtable(self, filename: str):
        """Save the Q-table to a file."""
        if self.q_table is None:
            raise ValueError("Q-table not built yet")
        
        np.savez(filename, 
                 q_table=self.q_table,
                 state_space_size=self.state_space_size,
                 action_space_size=self.action_space_size,
                 symbols=self.symbols,
                 symbol_probabilities=self.symbol_probabilities)
        print(f"Optimal Q-table saved to {filename}")
    
    def load_qtable(self, filename: str):
        """Load a Q-table from a file."""
        data = np.load(filename)
        self.q_table = data['q_table']
        print(f"Optimal Q-table loaded from {filename}")

class OptimalEVAgent:
    """Agent that uses the pre-computed optimal Q-table for decision making."""
    
    def __init__(self, game: SlotMachine, qtable: OptimalEVQTable):
        self.game = game
        self.qtable = qtable
        self.episode_rewards = []
        self.action_choices = {'cash_out': 0, 'respin': 0}
        
    def choose_action(self, state: Tuple[str, ...], respins_used: int) -> int:
        """Choose action using the optimal Q-table."""
        # Convert state to symbol counts
        symbol_counts = [state.count(symbol) for symbol in self.game.symbols]
        
        # Get optimal action
        action, value = self.qtable.get_optimal_action(symbol_counts, respins_used)
        
        # Track action choices
        if action == 8:
            self.action_choices['cash_out'] += 1
        else:
            self.action_choices['respin'] += 1
        
        return action
    
    def simulate_episode(self, max_steps: int = 20) -> float:
        """Simulate one episode using optimal strategy."""
        state = self.game.reset()
        total_reward = 0
        
        for step in range(max_steps):
            action = self.choose_action(self.game.current_state, self.game.respins_used)
            reward, next_state, done = self.game.step(action)
            total_reward += reward
            
            if done:
                break
        
        return total_reward
    
    def evaluate_strategy(self, num_episodes: int = 10000) -> dict:
        """Evaluate the optimal strategy."""
        print(f"Evaluating Optimal EV Strategy for {num_episodes} episodes...")
        
        rewards = []
        for episode in tqdm(range(num_episodes), desc="Evaluating Optimal Strategy"):
            reward = self.simulate_episode()
            rewards.append(reward)
            self.episode_rewards.append(reward)
        
        avg_reward = np.mean(rewards)
        std_reward = np.std(rewards)
        profitable_episodes = np.sum(np.array(rewards) > 0)
        
        return {
            'avg_reward': avg_reward,
            'std_reward': std_reward,
            'profitable_episodes': profitable_episodes,
            'profitable_rate': (profitable_episodes / num_episodes) * 100,
            'best_reward': np.max(rewards),
            'worst_reward': np.min(rewards),
            'cash_out_percentage': (self.action_choices['cash_out'] / sum(self.action_choices.values())) * 100,
            'respin_percentage': (self.action_choices['respin'] / sum(self.action_choices.values())) * 100
        }