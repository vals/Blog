class SlotMachine {
    constructor() {
        this.symbols = ['Blank', 'Coin', 'Coin Stack', 'Snake', 'Net', 'x2', 'Clover', 'Crown'];
        this.symbolProbabilities = [0.28, 0.30, 0.10, 0.10, 0.04, 0.09, 0.01, 0.08];
        this.numReels = 4;
        this.numSymbols = this.symbols.length;
        
        this.actions = {
            0: { type: 'respin_symbol', symbol: 'Blank', cost: 1, description: 'Respin a Blank' },
            1: { type: 'respin_symbol', symbol: 'Coin', cost: 1, description: 'Respin a Coin' },
            2: { type: 'respin_symbol', symbol: 'Coin Stack', cost: 1, description: 'Respin a Coin Stack' },
            3: { type: 'respin_symbol', symbol: 'Snake', cost: 1, description: 'Respin a Snake' },
            4: { type: 'respin_symbol', symbol: 'Net', cost: 1, description: 'Respin a Net' },
            5: { type: 'respin_symbol', symbol: 'x2', cost: 1, description: 'Respin a x2' },
            6: { type: 'respin_symbol', symbol: 'Clover', cost: 1, description: 'Respin a Clover' },
            7: { type: 'respin_symbol', symbol: 'Crown', cost: 1, description: 'Respin a Crown' },
            8: { type: 'cash_out', cost: 0, description: 'Cash out current winnings' }
        };
        
        this.maxRespins = 5;
        this.initialSpinCost = 1;
        
        this.currentState = this.generateRandomState();
        this.respinsUsed = 0;
        this.totalCost = 0;
        this.episodeActive = false;
        this.balance = 100;
        this.episodeStartBalance = 100;
    }
    
    generateRandomState() {
        const state = [];
        for (let i = 0; i < this.numReels; i++) {
            const random = Math.random();
            let cumulativeProb = 0;
            for (let j = 0; j < this.symbols.length; j++) {
                cumulativeProb += this.symbolProbabilities[j];
                if (random < cumulativeProb) {
                    state.push(this.symbols[j]);
                    break;
                }
            }
        }
        return state;
    }
    
    calculatePayout(state) {
        let payout = 0;
        
        // Count all symbols
        const symbolCounts = {};
        for (const symbol of this.symbols) {
            symbolCounts[symbol] = state.filter(s => s === symbol).length;
        }
        
        const coinCount = symbolCounts['Coin'] || 0;
        const coinStackCount = symbolCounts['Coin Stack'] || 0;
        const cloverCount = symbolCounts['Clover'] || 0;
        const netCount = symbolCounts['Net'] || 0;
        const crownCount = symbolCounts['Crown'] || 0;
        const snakeCount = symbolCounts['Snake'] || 0;
        
        // Check snake behavior - Net changes how Snake works
        if (snakeCount > 0 && netCount === 0) {
            // Snake with no Net: sets all winnings to zero
            return 0;
        }
        
        // Coin payouts
        if (coinCount === 3) {
            payout += 3;
        } else if (coinCount === 4) {
            payout += 5;
        }
        
        // Coin stack payouts (triple coins)
        if (coinStackCount === 3) {
            payout += 9;
        } else if (coinStackCount === 4) {
            payout += 15;
        }
        
        // Clover payouts - 10 gold each
        payout += cloverCount * 10;
        
        // Net payouts - 3 gold per snake when Net is present
        if (netCount > 0 && snakeCount > 0) {
            payout += netCount * snakeCount * 3;
        }
        
        // Crown jackpot
        if (crownCount === 4) {
            payout += 100;
        }
        
        // Apply 2x multiplier if present
        if (state.includes('x2')) {
            payout *= 2;
        }
        
        return payout;
    }
    
    step(action) {
        if (!(action in this.actions)) {
            throw new Error(`Invalid action: ${action}`);
        }
        
        const actionInfo = this.actions[action];
        
        if (actionInfo.type === 'cash_out') {
            // Cash out - collect current payout and end episode
            const payout = this.calculatePayout(this.currentState);
            const reward = payout - this.totalCost;
            this.balance += reward;
            this.episodeActive = false;
            
            return {
                reward: reward,
                state: [...this.currentState],
                done: true,
                payout: payout,
                totalCost: this.totalCost
            };
        } else if (actionInfo.type === 'respin_symbol') {
            // Check if respins are available
            if (this.respinsUsed >= this.maxRespins) {
                // No more respins allowed, must cash out
                const payout = this.calculatePayout(this.currentState);
                const reward = payout - this.totalCost;
                this.balance += reward;
                this.episodeActive = false;
                
                return {
                    reward: reward,
                    state: [...this.currentState],
                    done: true,
                    payout: payout,
                    totalCost: this.totalCost,
                    forcedCashOut: true
                };
            }
            
            // Check if action is valid (symbol exists)
            const symbolToRespin = actionInfo.symbol;
            if (!this.currentState.includes(symbolToRespin)) {
                // Invalid action - treat as no-op but still consume resources
                this.respinsUsed += 1;
                this.totalCost += actionInfo.cost;
                this.balance -= actionInfo.cost;
                
                return {
                    reward: -actionInfo.cost,
                    state: [...this.currentState],
                    done: false,
                    invalid: true
                };
            }
            
            // Find all positions with this symbol and randomly pick one
            const symbolPositions = [];
            for (let i = 0; i < this.currentState.length; i++) {
                if (this.currentState[i] === symbolToRespin) {
                    symbolPositions.push(i);
                }
            }
            
            if (symbolPositions.length === 0) {
                // This shouldn't happen due to the check above, but safety first
                return {
                    reward: 0,
                    state: [...this.currentState],
                    done: false,
                    invalid: true
                };
            }
            
            // Randomly select one position to respin
            const reelToRespin = symbolPositions[Math.floor(Math.random() * symbolPositions.length)];
            const cost = actionInfo.cost;
            
            // Respin the selected reel
            const newState = [...this.currentState];
            const random = Math.random();
            let cumulativeProb = 0;
            for (let j = 0; j < this.symbols.length; j++) {
                cumulativeProb += this.symbolProbabilities[j];
                if (random < cumulativeProb) {
                    newState[reelToRespin] = this.symbols[j];
                    break;
                }
            }
            
            this.currentState = newState;
            this.respinsUsed += 1;
            this.totalCost += cost;
            this.balance -= cost;
            
            // Episode continues, no immediate reward (reward comes at cash out)
            return {
                reward: 0,
                state: [...newState],
                done: false,
                cost: cost,
                reelRespun: reelToRespin
            };
        } else {
            throw new Error(`Unknown action type: ${actionInfo.type}`);
        }
    }
    
    // New method to respin a specific reel
    respinReel(reelIndex) {
        if (reelIndex < 0 || reelIndex >= this.numReels) {
            throw new Error(`Invalid reel index: ${reelIndex}`);
        }
        
        // Check if respins are available
        if (this.respinsUsed >= this.maxRespins) {
            // No more respins allowed, must cash out
            const payout = this.calculatePayout(this.currentState);
            const reward = payout - this.totalCost;
            this.balance += reward;
            this.episodeActive = false;
            
            return {
                reward: reward,
                state: [...this.currentState],
                done: true,
                payout: payout,
                totalCost: this.totalCost,
                forcedCashOut: true
            };
        }
        
        const cost = 1; // Respin cost
        const symbolBeforeRespin = this.currentState[reelIndex];
        
        // Respin the selected reel
        const newState = [...this.currentState];
        const random = Math.random();
        let cumulativeProb = 0;
        for (let j = 0; j < this.symbols.length; j++) {
            cumulativeProb += this.symbolProbabilities[j];
            if (random < cumulativeProb) {
                newState[reelIndex] = this.symbols[j];
                break;
            }
        }
        
        this.currentState = newState;
        this.respinsUsed += 1;
        this.totalCost += cost;
        this.balance -= cost;
        
        // Convert to equivalent action for data collection
        const equivalentAction = this.getActionForSymbol(symbolBeforeRespin);
        
        // Episode continues, no immediate reward (reward comes at cash out)
        return {
            reward: 0,
            state: [...newState],
            done: false,
            cost: cost,
            reelRespun: reelIndex,
            symbolRespun: symbolBeforeRespin,
            equivalentAction: equivalentAction
        };
    }
    
    // Helper method to get action ID for a symbol
    getActionForSymbol(symbol) {
        for (const [actionId, actionInfo] of Object.entries(this.actions)) {
            if (actionInfo.type === 'respin_symbol' && actionInfo.symbol === symbol) {
                return parseInt(actionId);
            }
        }
        return -1; // Not found
    }
    
    // New method to cash out
    cashOut() {
        return this.step(8); // Cash out action
    }
    
    reset() {
        // Start new episode with initial spin
        this.currentState = this.generateRandomState();
        this.respinsUsed = 0;
        this.totalCost = this.initialSpinCost; // Initial spin costs 1 Gold
        this.balance -= this.initialSpinCost;
        this.episodeActive = true;
        this.episodeStartBalance = this.balance + this.initialSpinCost;
        return [...this.currentState];
    }
    
    resetRun() {
        // Start new run with fresh balance
        this.balance = 100;
        this.reset();
        return [...this.currentState];
    }
    
    getValidActions() {
        const validActions = [];
        
        // Check if respins are still available
        if (this.respinsUsed >= this.maxRespins) {
            // Only cash out is available
            return [8]; // Cash out action
        }
        
        // Check each respin action
        for (const [actionId, actionInfo] of Object.entries(this.actions)) {
            const id = parseInt(actionId);
            if (actionInfo.type === 'cash_out') {
                validActions.push(id);
            } else if (actionInfo.type === 'respin_symbol') {
                const symbol = actionInfo.symbol;
                // Only valid if this symbol is present in current state
                if (this.currentState.includes(symbol)) {
                    validActions.push(id);
                }
            }
        }
        
        return validActions;
    }
    
    isActionValid(action) {
        return this.getValidActions().includes(action);
    }
    
    getStateAsCounts() {
        const counts = [];
        for (const symbol of this.symbols) {
            counts.push(this.currentState.filter(s => s === symbol).length);
        }
        counts.push(this.respinsUsed);
        return counts;
    }
    
    getCurrentPayout() {
        return this.calculatePayout(this.currentState);
    }
    
    getEpisodeBalance() {
        return this.balance - this.episodeStartBalance;
    }
    
    getGameState() {
        return {
            state: [...this.currentState],
            balance: this.balance,
            episodeBalance: this.getEpisodeBalance(),
            respinsUsed: this.respinsUsed,
            totalCost: this.totalCost,
            currentPayout: this.getCurrentPayout(),
            validActions: this.getValidActions(),
            episodeActive: this.episodeActive,
            stateCounts: this.getStateAsCounts()
        };
    }
}