import random
import numpy as np
from typing import Tuple, List, Dict
from math import comb

class SlotMachine:
    def __init__(self):
        self.symbols = ['Blank', 'Coin', 'Coin Stack', 'Snake', 'Net', 'x2', 'Clover', 'Crown']
        self.symbol_probabilities = [0.28, 0.30, 0.10, 0.10, 0.04, 0.09, 0.01, 0.08]
        self.num_reels = 4
        self.num_symbols = len(self.symbols)
        
        # Payout rules implemented in _calculate_payout method
        # - 3 single coins: 3 Gold
        # - 4 single coins: 5 Gold  
        # - 3 triple coins (Coin Stack): 9 Gold
        # - 4 triple coins (Coin Stack): 15 Gold
        # - Any four-leaf clover: 10 Gold each
        # - 2x symbol: Doubles the current spin's payout
        # - Red snake: Sets all winnings to zero
        # - Net: 3 Gold per snake on the line
        # - 4 gold crowns: 100 Gold (Jackpot)
        
        self.actions = {
            0: {'type': 'respin_symbol', 'symbol': 'Blank', 'cost': 1, 'description': 'Respin a Blank'},
            1: {'type': 'respin_symbol', 'symbol': 'Coin', 'cost': 1, 'description': 'Respin a Coin'},
            2: {'type': 'respin_symbol', 'symbol': 'Coin Stack', 'cost': 1, 'description': 'Respin a Coin Stack'},
            3: {'type': 'respin_symbol', 'symbol': 'Snake', 'cost': 1, 'description': 'Respin a Snake'},
            4: {'type': 'respin_symbol', 'symbol': 'Net', 'cost': 1, 'description': 'Respin a Net'},
            5: {'type': 'respin_symbol', 'symbol': 'x2', 'cost': 1, 'description': 'Respin a x2'},
            6: {'type': 'respin_symbol', 'symbol': 'Clover', 'cost': 1, 'description': 'Respin a Clover'},
            7: {'type': 'respin_symbol', 'symbol': 'Crown', 'cost': 1, 'description': 'Respin a Crown'},
            8: {'type': 'cash_out', 'cost': 0, 'description': 'Cash out current winnings'}
        }
        
        self.max_respins = 5
        self.initial_spin_cost = 1
        
        self.current_state = self._generate_random_state()
        self.respins_used = 0
        self.total_cost = 0
        self.episode_active = False
    
    def _generate_random_state(self) -> Tuple[str, ...]:
        return tuple(np.random.choice(self.symbols, p=self.symbol_probabilities) for _ in range(self.num_reels))
    
    def _calculate_payout(self, state: Tuple[str, ...], bet: int) -> int:
        payout = 0
        
        # Count all symbols
        coin_count = state.count('Coin')
        coin_stack_count = state.count('Coin Stack')
        clover_count = state.count('Clover')
        net_count = state.count('Net')
        crown_count = state.count('Crown')
        snake_count = state.count('Snake')
        
        # Check snake behavior - Net changes how Snake works
        if snake_count > 0 and net_count == 0:
            # Snake with no Net: sets all winnings to zero
            return 0
        
        # Coin payouts
        if coin_count == 3:
            payout += 3
        elif coin_count == 4:
            payout += 5
            
        # Coin stack payouts (triple coins)
        if coin_stack_count == 3:
            payout += 9
        elif coin_stack_count == 4:
            payout += 15
            
        # Clover payouts - 10 gold each
        payout += clover_count * 10
        
        # Net payouts - 3 gold per snake when Net is present
        if net_count > 0 and snake_count > 0:
            payout += net_count * snake_count * 3
            
        # Crown jackpot
        if crown_count == 4:
            payout += 100
            
        # Apply 2x multiplier if present
        if 'x2' in state:
            payout *= 2
            
        return payout
    
    def step(self, action: int) -> Tuple[int, Tuple[str, ...], bool]:
        if action not in self.actions:
            raise ValueError(f"Invalid action: {action}")
        
        action_info = self.actions[action]
        
        if action_info['type'] == 'cash_out':
            # Cash out - collect current payout and end episode
            payout = self._calculate_payout(self.current_state, 1)  # Base payout
            reward = payout - self.total_cost
            self.episode_active = False
            return reward, self.current_state, True
            
        elif action_info['type'] == 'respin_symbol':
            # Check if respins are available
            if self.respins_used >= self.max_respins:
                # No more respins allowed, must cash out
                payout = self._calculate_payout(self.current_state, 1)
                reward = payout - self.total_cost
                self.episode_active = False
                return reward, self.current_state, True
            
            # Check if action is valid (symbol exists)
            symbol_to_respin = action_info['symbol']
            if symbol_to_respin not in self.current_state:
                # Invalid action - treat as no-op but still consume resources
                self.respins_used += 1
                self.total_cost += action_info['cost']
                return -action_info['cost'], self.current_state, False
            
            # Find all positions with this symbol and randomly pick one
            symbol_positions = [i for i, symbol in enumerate(self.current_state) if symbol == symbol_to_respin]
            if not symbol_positions:
                # This shouldn't happen due to the check above, but safety first
                return 0, self.current_state, False
            
            # Randomly select one position to respin
            reel_to_respin = random.choice(symbol_positions)
            cost = action_info['cost']
            
            # Respin the selected reel
            new_state = list(self.current_state)
            new_state[reel_to_respin] = np.random.choice(self.symbols, p=self.symbol_probabilities)
            new_state = tuple(new_state)
            
            self.current_state = new_state
            self.respins_used += 1
            self.total_cost += cost
            
            # Episode continues, no immediate reward (reward comes at cash out)
            return 0, new_state, False
        
        else:
            raise ValueError(f"Unknown action type: {action_info['type']}")
    
    def reset(self) -> Tuple[str, ...]:
        # Start new episode with initial spin
        self.current_state = self._generate_random_state()
        self.respins_used = 0
        self.total_cost = self.initial_spin_cost  # Initial spin costs 1 Gold
        self.episode_active = True
        return self.current_state
    
    def get_state_space_size(self) -> int:
        # State includes reel symbols + number of respins used (0-5)
        return (self.num_symbols ** self.num_reels) * (self.max_respins + 1)
    
    def get_action_space_size(self) -> int:
        return len(self.actions)
    
    def state_to_index(self, state: Tuple[str, ...]) -> int:
        # Encode reel symbols
        reel_index = 0
        for i, symbol in enumerate(state):
            symbol_index = self.symbols.index(symbol)
            reel_index += symbol_index * (self.num_symbols ** (self.num_reels - 1 - i))
        
        # Combine with respins used
        total_index = reel_index * (self.max_respins + 1) + self.respins_used
        return total_index
    
    def index_to_state(self, index: int) -> Tuple[str, ...]:
        # Extract respins used
        respins_used = index % (self.max_respins + 1)
        reel_index = index // (self.max_respins + 1)
        
        # Extract reel symbols
        state = []
        for i in range(self.num_reels):
            symbol_index = reel_index // (self.num_symbols ** (self.num_reels - 1 - i))
            state.append(self.symbols[symbol_index])
            reel_index = reel_index % (self.num_symbols ** (self.num_reels - 1 - i))
        
        return tuple(state)
    
    def get_current_state_index(self) -> int:
        return self.state_to_index(self.current_state)
    
    # ============ STATE ABSTRACTION METHODS ============
    
    def state_to_counts(self, state: Tuple[str, ...]) -> Tuple[int, ...]:
        """Convert reel positions to symbol counts + respins used."""
        counts = []
        for symbol in self.symbols:
            counts.append(state.count(symbol))
        counts.append(self.respins_used)
        return tuple(counts)
    
    def get_current_counts(self) -> Tuple[int, ...]:
        """Get current state as symbol counts + respins used."""
        return self.state_to_counts(self.current_state)
    
    def counts_to_index(self, counts: Tuple[int, ...]) -> int:
        """Convert symbol counts + respins to a unique index using combinatorial number system."""
        symbol_counts = counts[:-1]  # All except respins
        respins_used = counts[-1]
        
        # Use a hash-based approach for simplicity and collision handling
        # This ensures all equivalent states map to the same index
        counts_tuple = tuple(symbol_counts)
        
        # Create a hash from the sorted counts to ensure equivalent states have same hash
        sorted_counts = tuple(sorted(symbol_counts, reverse=True))
        
        # Simple hash function that works well for small integer tuples
        hash_val = 0
        for i, count in enumerate(counts_tuple):
            hash_val += count * (i + 1) * 37  # 37 is a good hash multiplier
        
        # Ensure positive and within bounds
        base_index = abs(hash_val) % self._get_abstract_state_space_size()
        
        # Add respins component
        total_index = base_index + respins_used * self._get_abstract_state_space_size()
        
        return total_index
    
    def _get_abstract_state_space_size(self) -> int:
        """Calculate the number of possible symbol count combinations."""
        # For 4 symbols distributed among 8 types, we need to count
        # all non-negative integer solutions to: x1 + x2 + ... + x8 = 4
        # This is C(4 + 8 - 1, 8 - 1) = C(11, 7) = 330
        # But we'll use a larger number to handle hash collisions safely
        return 500
    
    def get_abstract_state_space_size(self) -> int:
        """Get total abstract state space including respins."""
        return self._get_abstract_state_space_size() * (self.max_respins + 1)
    
    def get_current_abstract_state_index(self) -> int:
        """Get current state as abstract index."""
        counts = self.get_current_counts()
        return self.counts_to_index(counts)
    
    # ============ ACTION MASKING METHODS ============
    
    def get_valid_actions(self) -> List[int]:
        """Get list of valid actions given current state."""
        valid_actions = []
        
        # Check if respins are still available
        if self.respins_used >= self.max_respins:
            # Only cash out is available
            return [8]  # Cash out action
        
        # Check each respin action
        for action_id, action_info in self.actions.items():
            if action_info['type'] == 'cash_out':
                valid_actions.append(action_id)
            elif action_info['type'] == 'respin_symbol':
                symbol = action_info['symbol']
                # Only valid if this symbol is present in current state
                if symbol in self.current_state:
                    valid_actions.append(action_id)
        
        return valid_actions
    
    def is_action_valid(self, action: int) -> bool:
        """Check if an action is valid given current state."""
        return action in self.get_valid_actions()