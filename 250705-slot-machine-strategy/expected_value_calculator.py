import numpy as np
from typing import Tuple, List, Dict
from math import comb
from game_rules import SlotMachine

class ExpectedValueCalculator:
    """Calculate exact expected values for slot machine actions based on game theory."""
    
    def __init__(self, game: SlotMachine):
        self.game = game
        self.symbols = game.symbols
        self.symbol_probabilities = game.symbol_probabilities
        self.symbol_to_prob = dict(zip(self.symbols, self.symbol_probabilities))
        
    def calculate_action_expected_value(self, state: Tuple[str, ...], action: int, respins_used: int) -> float:
        """Calculate expected value for a specific action in a given state."""
        if action == 8:  # Cash out
            return self._calculate_cash_out_value(state)
        else:
            return self._calculate_respin_expected_value(state, action, respins_used)
    
    def _calculate_cash_out_value(self, state: Tuple[str, ...]) -> float:
        """Calculate immediate cash-out value (no future uncertainty)."""
        return self.game._calculate_payout(state, 1)
    
    def _calculate_respin_expected_value(self, state: Tuple[str, ...], action: int, respins_used: int) -> float:
        """Calculate expected value of respinning a specific symbol type."""
        if respins_used >= self.game.max_respins:
            return float('-inf')  # No more respins allowed
            
        action_info = self.game.actions[action]
        symbol_to_respin = action_info['symbol']
        respin_cost = action_info['cost']
        
        # Check if symbol exists in current state
        if symbol_to_respin not in state:
            return float('-inf')  # Invalid action
        
        # Find positions of this symbol
        symbol_positions = [i for i, symbol in enumerate(state) if symbol == symbol_to_respin]
        
        # Calculate expected value for each possible position to respin
        total_expected_value = 0.0
        for position in symbol_positions:
            ev = self._calculate_position_respin_ev(state, position, respins_used)
            total_expected_value += ev
        
        # Average across all possible positions (since we randomly select one)
        average_ev = total_expected_value / len(symbol_positions)
        
        # Subtract the cost of respinning
        return average_ev - respin_cost
    
    def _calculate_position_respin_ev(self, state: Tuple[str, ...], position: int, respins_used: int) -> float:
        """Calculate expected value of respinning a specific position."""
        expected_value = 0.0
        
        # Consider all possible new symbols at this position
        for new_symbol, prob in zip(self.symbols, self.symbol_probabilities):
            new_state = list(state)
            new_state[position] = new_symbol
            new_state = tuple(new_state)
            
            # Calculate value of this new state
            state_value = self._calculate_state_value(new_state, respins_used + 1)
            expected_value += prob * state_value
            
        return expected_value
    
    def _calculate_state_value(self, state: Tuple[str, ...], respins_used: int) -> float:
        """Calculate the expected value of being in a specific state."""
        # Option 1: Cash out immediately
        cash_out_value = self._calculate_cash_out_value(state)
        
        # Option 2: Continue with optimal respin strategy (if respins available)
        if respins_used < self.game.max_respins:
            best_respin_value = float('-inf')
            
            # Check all possible respin actions
            for action in range(8):  # Actions 0-7 are respin actions
                action_info = self.game.actions[action]
                symbol_to_respin = action_info['symbol']
                
                if symbol_to_respin in state:
                    respin_ev = self._calculate_respin_expected_value(state, action, respins_used)
                    best_respin_value = max(best_respin_value, respin_ev)
            
            # If no valid respin actions, best_respin_value remains -inf
            if best_respin_value == float('-inf'):
                return cash_out_value
            else:
                return max(cash_out_value, best_respin_value)
        else:
            # No respins left, must cash out
            return cash_out_value
    
    def get_best_action(self, state: Tuple[str, ...], respins_used: int) -> Tuple[int, float]:
        """Get the action with highest expected value and its EV."""
        best_action = 8  # Default to cash out
        best_ev = self._calculate_cash_out_value(state)
        
        # Check all respin actions
        for action in range(8):
            action_info = self.game.actions[action]
            symbol_to_respin = action_info['symbol']
            
            if symbol_to_respin in state and respins_used < self.game.max_respins:
                ev = self._calculate_respin_expected_value(state, action, respins_used)
                if ev > best_ev:
                    best_ev = ev
                    best_action = action
        
        return best_action, best_ev
    
    def get_positive_ev_actions(self, state: Tuple[str, ...], respins_used: int) -> List[Tuple[int, float]]:
        """Get all actions with positive expected value."""
        cash_out_value = self._calculate_cash_out_value(state)
        positive_ev_actions = []
        
        # Check all respin actions
        for action in range(8):
            action_info = self.game.actions[action]
            symbol_to_respin = action_info['symbol']
            
            if symbol_to_respin in state and respins_used < self.game.max_respins:
                ev = self._calculate_respin_expected_value(state, action, respins_used)
                if ev > cash_out_value:  # Positive compared to cashing out
                    positive_ev_actions.append((action, ev))
        
        return positive_ev_actions
    
    def analyze_state_potential(self, state: Tuple[str, ...]) -> Dict[str, float]:
        """Analyze the potential of a state for different winning conditions."""
        counts = {}
        for symbol in self.symbols:
            counts[symbol] = state.count(symbol)
        
        analysis = {}
        
        # Jackpot potential (4 Crowns)
        crown_count = counts['Crown']
        if crown_count >= 1:
            # Probability of getting remaining crowns needed
            crowns_needed = 4 - crown_count
            prob_crown = self.symbol_to_prob['Crown']
            analysis['jackpot_potential'] = (prob_crown ** crowns_needed) * 100  # 100 Gold jackpot
        
        # Coin completion potential
        coin_count = counts['Coin']
        if coin_count >= 2:
            coins_needed = min(4, 4 - coin_count)  # Can get 3 or 4 coin payout
            prob_coin = self.symbol_to_prob['Coin']
            if coin_count == 2:
                analysis['coin_3_potential'] = prob_coin * 3  # 3 Gold for 3 coins
                analysis['coin_4_potential'] = (prob_coin ** 2) * 5  # 5 Gold for 4 coins
            elif coin_count == 3:
                analysis['coin_4_potential'] = prob_coin * 5  # 5 Gold for 4 coins
        
        # Coin stack potential
        stack_count = counts['Coin Stack']
        if stack_count >= 2:
            prob_stack = self.symbol_to_prob['Coin Stack']
            if stack_count == 2:
                analysis['stack_3_potential'] = prob_stack * 9  # 9 Gold for 3 stacks
                analysis['stack_4_potential'] = (prob_stack ** 2) * 15  # 15 Gold for 4 stacks
            elif stack_count == 3:
                analysis['stack_4_potential'] = prob_stack * 15  # 15 Gold for 4 stacks
        
        # Clover potential
        clover_count = counts['Clover']
        prob_clover = self.symbol_to_prob['Clover']
        analysis['clover_potential'] = (4 - clover_count) * prob_clover * 10  # 10 Gold per clover
        
        # Snake-Net interaction potential
        snake_count = counts['Snake']
        net_count = counts['Net']
        if snake_count > 0 and net_count > 0:
            # Current bonus from snake-net interaction
            analysis['snake_net_bonus'] = snake_count * net_count * 3
            # Potential for more snakes
            prob_snake = self.symbol_to_prob['Snake']
            analysis['additional_snake_potential'] = (4 - snake_count) * prob_snake * net_count * 3
        elif snake_count == 0 and net_count > 0:
            # Potential to get snakes when we have nets
            prob_snake = self.symbol_to_prob['Snake']
            analysis['snake_acquisition_potential'] = 4 * prob_snake * net_count * 3
        
        return analysis