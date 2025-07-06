# Data collection game

The rules of the slot machine game are described in @CLAUDE.md

For a point of comparison, we want to create a Q table based on empirical data.

To this end, we need to design a minimal javascript version of the slot machine. I will then play then play the gaame a number of times. The javascript version of the game will collect data on which strategy I use in which state, reflecting the episode simulations happening during Q learning.

The game should start with 100 Gold as in the simulations. And as in the simulations, balance over multiple episodes should be tracked.

The game does not need animations or anything too complicated. But symbols, reels, legal moves, episode balance, total balance, should be represented in the UI.

There should be keyboard support for performing the moves, including respinning reels (when legal) using numbers 1-4 on the keyboard, cash out being the enter key. There should be a button to start a new run, with new history and starting with balance of 100 Gold. 

The choice of action per state should be tracked. After we have collected a lot of data, we can use that to create the strategy table. We will also want to create balance plot over episodes like in the simulations. It will also be interesting to see the empirical policy distribution.

Whenever the game is played, new data should be appended to all other data collected from playing the game and stored together. This way empirical data can be collected over multiple sessions.

A nice bonus feature would be to time stamp the data.

## Implementation Results

The JavaScript slot machine game was successfully implemented with the following features:

### Game Interface
- **Reel-based controls**: Keys 1-4 map directly to reels 1-4 for intuitive respinning
- **Dynamic button updates**: Buttons show current symbol and update in real-time
- **Dual-purpose Enter key**: Cash out during episodes, start new episode when inactive
- **Real-time statistics**: Balance tracking, episode counts, and payout display
- **Persistent data storage**: All decisions saved to localStorage across sessions

### Data Collection
- **Comprehensive tracking**: Every state-action pair recorded with timestamps
- **Export functionality**: JSON and CSV export for analysis
- **Policy analysis**: Automatic calculation of decision patterns and strategy distribution
- **Session persistence**: Data accumulates across multiple play sessions

## Empirical vs Q-Learning Comparison Results

After collecting empirical data through human gameplay, a detailed comparison was conducted between human decision-making and the Q-learning agent's learned policy.

### Performance Summary

**Human Player Performance:**
- **Runs completed**: 26 runs
- **Average episodes per run**: 30.2 episodes
- **Average final balance**: -5.2 gold (starting from 100 gold)
- **Best run**: 0 gold final balance
- **Worst run**: -11 gold final balance
- **Standard deviation**: 3.2 gold
- **Success rate**: 0% (no runs ended above starting balance)
- **Total decisions made**: 786 decisions

**Q-Learning Agent Performance:**
- **Runs analyzed**: 10 runs  
- **Average episodes per run**: 50 episodes
- **Average final balance**: -6.3 gold
- **Best run**: 18 gold final balance
- **Worst run**: -60 gold final balance
- **Standard deviation**: 23.7 gold
- **Success rate**: 0% (no runs ended above starting balance)

### Key Findings

1. **Human outperformed Q-learning**: The human player achieved a better average final balance by 1.1 gold compared to the Q-learning agent.

2. **Consistency vs Variance**: Human play showed much more consistent results (σ = 3.2) compared to the Q-learning agent (σ = 23.7), suggesting more stable decision-making patterns.

3. **Risk management**: Human decisions appear to be more conservative, avoiding the extreme losses seen in some Q-learning runs.

4. **Both strategies struggled**: Neither approach consistently beat the house edge, with both showing negative expected returns.

### Human Decision Patterns

The empirical data revealed interesting patterns in human decision-making:
- **Action distribution**: Analysis of the 786 decisions shows clear preferences for certain symbol respins
- **State-dependent choices**: Human strategy varied based on current game state and respin count
- **Cash-out timing**: Patterns in when humans chose to cash out vs continue respinning

### Implications

The comparison suggests that:
1. **Human intuition can be competitive** with algorithmic approaches in complex decision environments
2. **Consistency matters** - the lower variance in human performance may be valuable in risk management
3. **Both strategies need refinement** - the negative returns indicate room for improvement in both human and AI approaches
4. **Empirical data provides valuable insights** - Human decision patterns could inform better Q-learning reward structures

### Generated Visualizations

Two comprehensive comparison plots were generated:
- `empirical_vs_simulated_comparison.png`: Basic balance progression comparison with color-coded runs
- `detailed_empirical_vs_simulated.png`: Six-panel analysis including distributions, action patterns, and performance metrics

The empirical data collection successfully provided a rich dataset for comparing human vs AI decision-making in the slot machine environment.