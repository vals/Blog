class GameUI {
    constructor(game, dataCollector) {
        this.game = game;
        this.dataCollector = dataCollector;
        this.logMessages = [];
        this.maxLogMessages = 100;
        
        this.initializeUI();
        this.setupEventListeners();
        this.updateDisplay();
    }
    
    initializeUI() {
        // Get all UI elements
        this.elements = {
            reels: [
                document.getElementById('reelDisplay0'),
                document.getElementById('reelDisplay1'),
                document.getElementById('reelDisplay2'),
                document.getElementById('reelDisplay3')
            ],
            reelSymbols: [
                document.querySelector('#reelDisplay0 .reel-symbol'),
                document.querySelector('#reelDisplay1 .reel-symbol'),
                document.querySelector('#reelDisplay2 .reel-symbol'),
                document.querySelector('#reelDisplay3 .reel-symbol')
            ],
            currentBalance: document.getElementById('currentBalance'),
            episodeBalance: document.getElementById('episodeBalance'),
            totalBalance: document.getElementById('totalBalance'),
            respinsUsed: document.getElementById('respinsUsed'),
            episodeCount: document.getElementById('episodeCount'),
            currentPayout: document.getElementById('currentPayout'),
            gameLog: document.getElementById('gameLog'),
            
            // Reel buttons
            reelButtons: [
                document.getElementById('reel0'),
                document.getElementById('reel1'),
                document.getElementById('reel2'),
                document.getElementById('reel3')
            ],
            cashOutBtn: document.getElementById('cashOut'),
            
            // Control buttons
            newEpisodeBtn: document.getElementById('newEpisodeBtn'),
            newRunBtn: document.getElementById('newRunBtn'),
            exportDataBtn: document.getElementById('exportDataBtn'),
            clearDataBtn: document.getElementById('clearDataBtn'),
            
            // Data export buttons
            exportJsonBtn: document.getElementById('exportJsonBtn'),
            exportCsvBtn: document.getElementById('exportCsvBtn'),
            showStatsBtn: document.getElementById('showStatsBtn'),
            
            // Data stats
            totalDecisions: document.getElementById('totalDecisions'),
            totalEpisodes: document.getElementById('totalEpisodes'),
            totalRuns: document.getElementById('totalRuns'),
            avgReward: document.getElementById('avgReward')
        };
    }
    
    setupEventListeners() {
        // Reel button listeners
        this.elements.reelButtons.forEach((button, index) => {
            if (button) {
                button.addEventListener('click', () => this.handleReelAction(index));
            }
        });
        
        // Cash out button listener
        if (this.elements.cashOutBtn) {
            this.elements.cashOutBtn.addEventListener('click', () => {
                if (this.game.episodeActive) {
                    this.handleCashOut();
                } else {
                    this.startNewEpisode();
                }
            });
        }
        
        // Control button listeners
        this.elements.newEpisodeBtn.addEventListener('click', () => this.startNewEpisode());
        this.elements.newRunBtn.addEventListener('click', () => this.startNewRun());
        this.elements.exportDataBtn.addEventListener('click', () => this.exportData());
        this.elements.clearDataBtn.addEventListener('click', () => this.clearData());
        
        // Data export button listeners
        this.elements.exportJsonBtn.addEventListener('click', () => this.exportJSON());
        this.elements.exportCsvBtn.addEventListener('click', () => this.exportCSV());
        this.elements.showStatsBtn.addEventListener('click', () => this.showStatistics());
        
        // Keyboard listeners
        document.addEventListener('keydown', (event) => this.handleKeyboard(event));
    }
    
    handleKeyboard(event) {
        // Prevent default behavior for game keys
        const gameKeys = ['1', '2', '3', '4', 'Enter', ' ', 'r', 'R'];
        if (gameKeys.includes(event.key)) {
            event.preventDefault();
        }
        
        switch (event.key) {
            case '1':
                this.handleReelAction(0);
                break;
            case '2':
                this.handleReelAction(1);
                break;
            case '3':
                this.handleReelAction(2);
                break;
            case '4':
                this.handleReelAction(3);
                break;
            case 'Enter':
                if (this.game.episodeActive) {
                    this.handleCashOut();
                } else {
                    this.startNewEpisode();
                }
                break;
            case ' ':
                this.startNewEpisode();
                break;
            case 'r':
            case 'R':
                this.startNewRun();
                break;
        }
    }
    
    handleReelAction(reelIndex) {
        if (!this.game.episodeActive) {
            this.addLogMessage('Episode not active. Start a new episode.', 'error');
            return;
        }
        
        const gameState = this.game.getGameState();
        const symbolAtReel = gameState.state[reelIndex];
        
        // Log the reel action
        this.addLogMessage(`Respin Reel ${reelIndex + 1} (${symbolAtReel})`, 'action');
        
        // Execute the reel respin
        const result = this.game.respinReel(reelIndex);
        
        if (result.done) {
            // Episode finished due to respin limit
            this.dataCollector.finishCurrentEpisode(this.game.balance, result.reward);
            
            const message = result.forcedCashOut ? 
                `Forced cash out! Reward: ${result.reward} Gold` :
                `Episode ended. Reward: ${result.reward} Gold`;
            this.addLogMessage(message, 'reward');
            
            this.addLogMessage(`Episode ${this.dataCollector.episodeCounter} finished.`, 'episode');
        } else if (result.cost) {
            // Record the decision with equivalent action for data compatibility
            this.dataCollector.recordDecision(gameState.state, result.equivalentAction, gameState);
            this.addLogMessage(`Respin cost: ${result.cost} Gold`, 'action');
        }
        
        this.updateDisplay();
    }
    
    handleCashOut() {
        if (!this.game.episodeActive) {
            this.addLogMessage('Episode not active. Start a new episode.', 'error');
            return;
        }
        
        const gameState = this.game.getGameState();
        
        // Record the cash out decision
        this.dataCollector.recordDecision(gameState.state, 8, gameState);
        
        // Execute the cash out
        const result = this.game.cashOut();
        
        // Log the action
        this.addLogMessage(`Cashed out! Reward: ${result.reward} Gold`, 'reward');
        
        // Episode finished
        this.dataCollector.finishCurrentEpisode(this.game.balance, result.reward);
        this.addLogMessage(`Episode ${this.dataCollector.episodeCounter} finished.`, 'episode');
        
        this.updateDisplay();
    }
    
    startNewEpisode() {
        if (this.game.balance < 1) {
            this.addLogMessage('Insufficient balance for new episode!', 'error');
            return;
        }
        
        this.dataCollector.startNewEpisode(this.game.balance);
        this.game.reset();
        this.addLogMessage(`New episode started. Episode ${this.dataCollector.episodeCounter}`, 'episode');
        this.updateDisplay();
    }
    
    startNewRun() {
        if (this.game.episodeActive) {
            const confirmed = confirm('End current episode and start new run?');
            if (!confirmed) return;
        }
        
        // Finish current episode if active
        if (this.game.episodeActive) {
            const gameState = this.game.getGameState();
            this.dataCollector.finishCurrentEpisode(this.game.balance, gameState.episodeBalance);
        }
        
        // Finish current run
        this.dataCollector.finishCurrentRun(this.game.balance);
        
        // Start new run
        this.dataCollector.startNewRun();
        this.game.resetRun();
        this.addLogMessage(`New run started. Run ${this.dataCollector.runCounter}`, 'episode');
        
        // Start first episode of new run
        this.dataCollector.startNewEpisode(this.game.balance);
        this.game.reset();
        this.addLogMessage(`New episode started. Episode ${this.dataCollector.episodeCounter}`, 'episode');
        
        this.updateDisplay();
    }
    
    exportData() {
        const data = this.dataCollector.exportToJSON();
        this.downloadFile(data, 'slot_machine_data.json', 'application/json');
        this.addLogMessage('Data exported to JSON file', 'info');
    }
    
    exportJSON() {
        const data = this.dataCollector.exportToJSON();
        this.downloadFile(data, 'slot_machine_data.json', 'application/json');
        this.addLogMessage('Data exported to JSON file', 'info');
    }
    
    exportCSV() {
        const data = this.dataCollector.exportToCSV();
        this.downloadFile(data, 'slot_machine_data.csv', 'text/csv');
        this.addLogMessage('Data exported to CSV file', 'info');
    }
    
    showStatistics() {
        const stats = this.dataCollector.getStatistics();
        const policy = this.dataCollector.getPolicyDistribution();
        
        let message = 'STATISTICS:\\n';
        message += `Total Decisions: ${stats.totalDecisions}\\n`;
        message += `Total Episodes: ${stats.totalEpisodes}\\n`;
        message += `Total Runs: ${stats.totalRuns}\\n`;
        message += `Avg Episode Reward: ${stats.avgEpisodeReward.toFixed(2)}\\n`;
        message += `Avg Run Reward: ${stats.avgRunReward.toFixed(2)}\\n`;
        message += `Avg Decisions/Episode: ${stats.avgDecisionsPerEpisode.toFixed(2)}\\n`;
        message += `Avg Episodes/Run: ${stats.avgEpisodesPerRun.toFixed(2)}\\n\\n`;
        
        message += 'ACTION DISTRIBUTION:\\n';
        for (const [action, count] of Object.entries(stats.actionDistribution)) {
            const percentage = (count / stats.totalDecisions * 100).toFixed(1);
            message += `${this.dataCollector.getActionDescription(action)}: ${count} (${percentage}%)\\n`;
        }
        
        message += `\\nUNIQUE STATES VISITED: ${policy.length}\\n`;
        message += `SESSION DURATION: ${Math.round(stats.sessionDuration / 60000)} minutes`;
        
        alert(message);
    }
    
    clearData() {
        const confirmed = confirm('Are you sure you want to clear all data? This cannot be undone.');
        if (confirmed) {
            this.dataCollector.clearData();
            this.addLogMessage('All data cleared', 'info');
            this.updateDisplay();
        }
    }
    
    downloadFile(content, filename, contentType) {
        const blob = new Blob([content], { type: contentType });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = filename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }
    
    updateDisplay() {
        const gameState = this.game.getGameState();
        const stats = this.dataCollector.getStatistics();
        
        // Update reel symbols
        this.elements.reelSymbols.forEach((symbolElement, index) => {
            if (symbolElement) {
                symbolElement.textContent = gameState.state[index];
            }
        });
        
        // Update reel displays
        this.elements.reels.forEach((reel, index) => {
            if (reel) {
                reel.className = 'reel';
            }
        });
        
        // Update balance stats
        this.elements.currentBalance.textContent = gameState.balance;
        this.elements.episodeBalance.textContent = gameState.episodeBalance;
        this.elements.totalBalance.textContent = gameState.balance;
        this.elements.respinsUsed.textContent = gameState.respinsUsed;
        this.elements.episodeCount.textContent = this.dataCollector.episodeCounter;
        this.elements.currentPayout.textContent = gameState.currentPayout;
        
        // Update reel buttons - only enable if episode is active
        this.elements.reelButtons.forEach((button, index) => {
            if (button) {
                button.disabled = !gameState.episodeActive;
                button.className = gameState.episodeActive ? 'action-btn' : 'action-btn';
                
                // Update button text to show current symbol
                const symbol = gameState.state[index];
                button.textContent = `${index + 1}: Respin ${symbol}`;
            }
        });
        
        // Update cash out button
        if (this.elements.cashOutBtn) {
            // Button is always enabled - either for cash out or new episode
            this.elements.cashOutBtn.disabled = false;
            this.elements.cashOutBtn.className = 'action-btn cash-out';
            
            // Update button text based on episode state
            if (gameState.episodeActive) {
                this.elements.cashOutBtn.textContent = 'ENTER: Cash Out';
            } else {
                this.elements.cashOutBtn.textContent = 'ENTER: New Episode';
            }
        }
        
        // Update control buttons
        this.elements.newEpisodeBtn.disabled = gameState.episodeActive || gameState.balance < 1;
        
        // Also disable cash out/new episode button if insufficient balance and no active episode
        if (this.elements.cashOutBtn && !gameState.episodeActive && gameState.balance < 1) {
            this.elements.cashOutBtn.disabled = true;
            this.elements.cashOutBtn.textContent = 'ENTER: Insufficient Balance';
        }
        
        // Update data statistics
        this.elements.totalDecisions.textContent = stats.totalDecisions;
        this.elements.totalEpisodes.textContent = stats.totalEpisodes;
        this.elements.totalRuns.textContent = stats.totalRuns;
        this.elements.avgReward.textContent = stats.avgEpisodeReward.toFixed(2);
    }
    
    addLogMessage(message, type = 'info') {
        const timestamp = new Date().toLocaleTimeString();
        this.logMessages.push({ message, type, timestamp });
        
        // Keep only recent messages
        if (this.logMessages.length > this.maxLogMessages) {
            this.logMessages.shift();
        }
        
        this.updateGameLog();
    }
    
    updateGameLog() {
        const logElement = this.elements.gameLog;
        logElement.innerHTML = '';
        
        this.logMessages.forEach(log => {
            const entry = document.createElement('div');
            entry.className = `log-entry ${log.type}`;
            entry.textContent = `[${log.timestamp}] ${log.message}`;
            logElement.appendChild(entry);
        });
        
        // Scroll to bottom
        logElement.scrollTop = logElement.scrollHeight;
    }
}