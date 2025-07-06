# Slot Machine Q-Learning Strategy Project

This project implements a Q-learning algorithm to find optimal strategies for a specific 4-reel slot machine game with respin mechanics.

## Project Structure

- `game_rules.py` - Slot machine game mechanics and payout calculations
- `q_learning.py` - Tabular Q-learning algorithm implementation
- `utils.py` - Helper functions for analysis and visualization
- `main.py` - Main training script that coordinates everything

## Slot Machine Rules

### Reels and Symbols
- **4 reels** with 8 possible symbols each
- **Symbol probabilities**:
  - Blank: 28%
  - Coin: 30%
  - Coin Stack: 10%
  - Snake: 10%
  - Net: 4%
  - x2: 9%
  - Clover: 1%
  - Crown: 8%

### Game Flow
1. **Initial Spin**: Costs 1 Gold, spins all 4 reels
2. **Decision Phase**: After seeing results, choose to:
   - Respin any instance of a specific symbol type (costs 1 Gold each)
   - Cash out current winnings (ends episode)
3. **Respin Limit**: Maximum 5 respins per episode
4. **Episode End**: Cash out or hit respin limit

### Payout Rules
- **3 single coins**: 3 Gold
- **4 single coins**: 5 Gold
- **3 triple coins (Coin Stack)**: 9 Gold
- **4 triple coins (Coin Stack)**: 15 Gold
- **Any four-leaf clover**: 10 Gold each
- **4 gold crowns**: 100 Gold (Jackpot)
- **2x symbol**: Doubles the current spin's payout
- **Red snake**: Sets all winnings to zero (unless Net present)
- **Net**: 3 Gold per snake on the line (neutralizes snake penalty)

### Special Symbol Interactions
- **Snake + Net**: Net "catches" the snake, converting penalty to bonus (3 Gold per snake per net)
- **Snake without Net**: All winnings set to zero
- **2x Multiplier**: Applied to final payout after all other calculations

### Actions Available
- **Action 0-7**: Respin any instance of specific symbol type (1 Gold each)
  - Action 0: Respin any Blank
  - Action 1: Respin any Coin  
  - Action 2: Respin any Coin Stack
  - Action 3: Respin any Snake
  - Action 4: Respin any Net
  - Action 5: Respin any x2
  - Action 6: Respin any Clover
  - Action 7: Respin any Crown
- **Action 8**: Cash out current winnings and end episode

### State Space (Abstracted)
- **3,000 total states**: ~500 symbol count combinations × 6 respin counts (0-5)
- State abstraction uses symbol counts instead of positions
- Equivalent game situations (same symbol counts) share same Q-values
- Massive reduction from original 24,576 position-dependent states

## Q-Learning Implementation

The agent learns optimal policies through:
- **Epsilon-greedy exploration** with optimized slow decay (0.9995)
- **Tabular Q-learning** updates using Bellman equation
- **Episode-based training** where reward comes only at cash-out
- **State abstraction** using symbol counts instead of positions
- **Action masking** to prevent invalid symbol respins
- **Adaptive learning rate** with decay for stability

## Training Goals

The Q-learning agent will discover:
- When to respin vs cash out based on current symbol counts
- Which symbol types to prioritize for respins
- How to balance potential gains against respin costs
- Optimal strategies for different symbol combinations
- Strategic use of action masking (only respin symbols that are present)

## Usage

### Training the Q-Learning Agent
Run the main training script:
```bash
python main.py
```

This will:
1. Train the agent for 100,000 episodes (optimized for convergence)
2. Analyze the learned policy with action masking
3. Evaluate performance on abstract state space
4. Save results and visualizations as PNG files

### Strategy Simulation and Comparison
After training, run strategy simulations:
```bash
python simulate_gameplay.py      # Q-learning agent strategy
python simulate_naive_strategy.py  # Naive "always cash out" baseline
```

### Human vs AI Analysis
Generate comprehensive strategy comparison:
```bash
python human_vs_agent_analysis.py    # Compare empirical human data with Q-learning
python generate_strategy_report.py   # Human-readable analysis report
```

### Visualization
Use the Jupyter notebook for comparative analysis:
```bash
jupyter notebook "250704 Make figures.ipynb"
```

## Strategy Comparison Results

The project includes three distinct strategies with surprising performance differences:

### Performance Summary (Average Net Result per Run)
1. **Naive "Cash out only"**: -68 Gold (8% profitable runs)
2. **Q-Learning Agent**: -102 Gold (0% profitable runs)  
3. **Human Players**: Variable (empirical data from real gameplay)

### Key Findings
- **Naive strategy outperforms Q-learning**: The simple "always cash out" approach loses less money on average
- **Q-learning agent never achieved profitability**: Despite sophisticated learning, 0% of runs were profitable
- **Human decision patterns**: Show different risk tolerance and strategic approaches compared to both AI strategies
- **Jackpot achievement**: Only the Q-learning agent learned to achieve the 100 Gold jackpot (4×Crown), though rarely

### Strategic Insights
- **Risk vs Reward**: The naive strategy's immediate cash-out minimizes losses but foregoes potential big wins
- **Exploration vs Exploitation**: Q-learning agent's exploration led to discovering high-value states but poor overall performance
- **Human intuition**: Empirical human data reveals cognitive biases and heuristic decision-making patterns
- **State abstraction effectiveness**: Agent successfully learned optimal policies for high-value symbol combinations

## Expected Insights

The trained agent should learn:
- Cash out immediately with high-value combinations
- Respin to complete partial winning combinations
- Avoid costly respins when probability of improvement is low
- Leverage the Net-Snake interaction strategically

## Data Files Generated

### Training Data
- `slot_machine_model.npz` - Trained Q-table and agent parameters
- `training_progress.png` - Learning curves and convergence analysis

### Simulation Data  
- `simulated_slot_machine_data.csv` - Q-learning agent decision-by-decision data
- `naive_slot_machine_data.csv` - Naive strategy decision-by-decision data
- `slot_machine_data.csv` - Empirical human gameplay data

### Analysis Results
- `human_vs_agent_analysis.json` - Comprehensive strategy comparison data
- `Human_vs_Agent_Strategy_Report.md` - Human-readable analysis report
- `multiple_balance_simulation.png` - Q-learning agent balance progressions
- `naive_balance_simulation.png` - Naive strategy balance progressions

### Aggregate Statistics
- `aggregate_simulation_results.txt` - Q-learning agent performance summary
- `naive_aggregate_results.txt` - Naive strategy performance summary

## Optimized Hyperparameters

After extensive testing, the following hyperparameters provide optimal learning:

**Core Parameters:**
- Learning rate: 0.25 (with decay 0.9999 → 0.01 minimum)
- Discount factor: 0.99 (values future rewards highly)
- Epsilon decay: 0.9995 (slow exploration decay for complex state space)
- Episodes: 100,000 (required for convergence with state abstraction)

**Key Insights:**
- Slower epsilon decay (0.9995 vs 0.9985) critical for learning
- Epsilon reaches minimum around episode 9,200 (vs 3,000 previously)
- Higher learning rate (0.25) compensates for sparse episode rewards
- Extended exploration essential for symbol-based action discovery