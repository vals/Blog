#!/usr/bin/env python3

"""
Test script to validate the refactored slot machine system with state abstraction and symbol-based actions.
"""

from game_rules import SlotMachine
from q_learning import QLearningAgent

def test_state_abstraction():
    """Test that equivalent states have the same abstract representation."""
    print("Testing state abstraction...")
    
    game = SlotMachine()
    
    # Create two states with same symbol counts but different positions
    state1 = ('Coin', 'Blank', 'Coin', 'Blank')
    state2 = ('Blank', 'Coin', 'Blank', 'Coin')
    
    # Set the states and check their abstract representations
    game.current_state = state1
    game.respins_used = 0
    abstract_index1 = game.get_current_abstract_state_index()
    
    game.current_state = state2
    game.respins_used = 0
    abstract_index2 = game.get_current_abstract_state_index()
    
    print(f"State 1: {state1} -> Abstract index: {abstract_index1}")
    print(f"State 2: {state2} -> Abstract index: {abstract_index2}")
    print(f"Same abstract index: {abstract_index1 == abstract_index2}")
    
    # Test symbol counts
    counts1 = game.state_to_counts(state1)
    counts2 = game.state_to_counts(state2)
    print(f"Counts 1: {counts1}")
    print(f"Counts 2: {counts2}")
    print(f"Same counts: {counts1[:-1] == counts2[:-1]}")  # Exclude respins
    
    return abstract_index1 == abstract_index2

def test_action_masking():
    """Test that action masking works correctly."""
    print("\nTesting action masking...")
    
    game = SlotMachine()
    
    # Test state with some symbols
    test_state = ('Coin', 'Snake', 'Blank', 'Crown')
    game.current_state = test_state
    game.respins_used = 0
    
    valid_actions = game.get_valid_actions()
    print(f"Test state: {test_state}")
    print(f"Valid actions: {valid_actions}")
    
    # Check that we can respin symbols that are present
    expected_symbols = ['Coin', 'Snake', 'Blank', 'Crown']
    for action_id, action_info in game.actions.items():
        if action_info['type'] == 'respin_symbol':
            symbol = action_info['symbol']
            should_be_valid = symbol in expected_symbols
            is_valid = action_id in valid_actions
            print(f"  Action {action_id} (respin {symbol}): Expected {should_be_valid}, Got {is_valid}")
            if should_be_valid != is_valid:
                return False
    
    # Test with max respins used
    game.respins_used = game.max_respins
    valid_actions_max = game.get_valid_actions()
    print(f"With max respins used, valid actions: {valid_actions_max}")
    print(f"Only cash out available: {valid_actions_max == [8]}")
    
    return True

def test_symbol_based_respins():
    """Test that symbol-based respins work correctly."""
    print("\nTesting symbol-based respins...")
    
    game = SlotMachine()
    
    # Set a test state
    initial_state = ('Snake', 'Snake', 'Coin', 'Blank')
    game.current_state = initial_state
    game.respins_used = 0
    game.total_cost = 1
    game.episode_active = True
    
    print(f"Initial state: {initial_state}")
    
    # Try to respin a Snake (action 3)
    reward, new_state, done = game.step(3)  # Respin Snake
    
    print(f"After respinning Snake:")
    print(f"  New state: {new_state}")
    print(f"  Reward: {reward}")
    print(f"  Done: {done}")
    print(f"  Respins used: {game.respins_used}")
    print(f"  Total cost: {game.total_cost}")
    
    # Check that exactly one Snake was replaced
    snake_count_before = initial_state.count('Snake')
    snake_count_after = new_state.count('Snake')
    print(f"Snake count before: {snake_count_before}, after: {snake_count_after}")
    print(f"Exactly one Snake removed: {snake_count_before - snake_count_after == 1}")
    
    return snake_count_before - snake_count_after <= 1  # Could be 0 if respun to Snake again

def test_basic_training():
    """Test that the Q-learning agent can train with the new system."""
    print("\nTesting basic Q-learning training...")
    
    game = SlotMachine()
    agent = QLearningAgent(
        state_space_size=game.get_abstract_state_space_size(),
        action_space_size=game.get_action_space_size(),
        learning_rate=0.1,
        discount_factor=0.9,
        epsilon=0.5,
        epsilon_decay=0.99,
        epsilon_min=0.1
    )
    
    print(f"State space size: {game.get_abstract_state_space_size()}")
    print(f"Action space size: {game.get_action_space_size()}")
    
    # Run a few training episodes
    initial_epsilon = agent.epsilon
    rewards = agent.train(game, num_episodes=100, verbose=False, report_interval=50)
    final_epsilon = agent.epsilon
    
    print(f"Training completed:")
    print(f"  Episodes trained: {len(rewards)}")
    print(f"  Average reward: {sum(rewards)/len(rewards):.3f}")
    print(f"  Epsilon decay: {initial_epsilon:.3f} -> {final_epsilon:.3f}")
    
    return len(rewards) == 100

def main():
    """Run all tests."""
    print("Running refactored system validation tests...")
    print("=" * 50)
    
    tests = [
        test_state_abstraction,
        test_action_masking,
        test_symbol_based_respins,
        test_basic_training
    ]
    
    results = []
    for test in tests:
        try:
            result = test()
            results.append(result)
            print(f"✓ {test.__name__}: {'PASSED' if result else 'FAILED'}")
        except Exception as e:
            results.append(False)
            print(f"✗ {test.__name__}: ERROR - {e}")
    
    print("\n" + "=" * 50)
    print(f"Tests passed: {sum(results)}/{len(results)}")
    
    if all(results):
        print("🎉 All tests passed! The refactored system is working correctly.")
    else:
        print("⚠️  Some tests failed. Please check the implementation.")

if __name__ == "__main__":
    main()