class DataCollector {
    constructor() {
        this.data = {
            decisions: [],
            episodes: [],
            runs: [],
            sessionStart: Date.now()
        };
        this.currentEpisode = {
            decisions: [],
            startTime: null,
            endTime: null,
            startBalance: 0,
            endBalance: 0,
            reward: 0,
            runId: null,
            episodeId: null
        };
        this.currentRun = {
            episodes: [],
            startTime: null,
            endTime: null,
            startBalance: 100,
            endBalance: 0,
            runId: null
        };
        this.runCounter = 0;
        this.episodeCounter = 0;
        this.decisionCounter = 0;
        
        // Load existing data from localStorage
        this.loadData();
    }
    
    startNewRun() {
        // Save current run if it exists
        if (this.currentRun.runId !== null) {
            this.finishCurrentRun();
        }
        
        this.runCounter++;
        this.currentRun = {
            episodes: [],
            startTime: Date.now(),
            endTime: null,
            startBalance: 100,
            endBalance: 0,
            runId: this.runCounter
        };
        
        this.episodeCounter = 0; // Reset episode counter for new run
    }
    
    startNewEpisode(startBalance) {
        // Save current episode if it exists
        if (this.currentEpisode.episodeId !== null) {
            this.finishCurrentEpisode();
        }
        
        this.episodeCounter++;
        this.currentEpisode = {
            decisions: [],
            startTime: Date.now(),
            endTime: null,
            startBalance: startBalance,
            endBalance: 0,
            reward: 0,
            runId: this.currentRun.runId,
            episodeId: this.episodeCounter
        };
    }
    
    recordDecision(state, action, gameState, timestamp = null) {
        if (timestamp === null) {
            timestamp = Date.now();
        }
        
        this.decisionCounter++;
        const decision = {
            decisionId: this.decisionCounter,
            runId: this.currentRun.runId,
            episodeId: this.currentEpisode.episodeId,
            timestamp: timestamp,
            state: [...state],
            stateCounts: [...gameState.stateCounts],
            action: action,
            actionDescription: this.getActionDescription(action),
            respinsUsed: gameState.respinsUsed,
            balance: gameState.balance,
            currentPayout: gameState.currentPayout,
            totalCost: gameState.totalCost,
            validActions: [...gameState.validActions]
        };
        
        this.currentEpisode.decisions.push(decision);
        this.data.decisions.push(decision);
        
        // Save to localStorage after each decision
        this.saveData();
        
        return decision;
    }
    
    finishCurrentEpisode(endBalance = 0, reward = 0) {
        if (this.currentEpisode.episodeId === null) {
            return;
        }
        
        this.currentEpisode.endTime = Date.now();
        this.currentEpisode.endBalance = endBalance;
        this.currentEpisode.reward = reward;
        
        const episodeData = {
            ...this.currentEpisode,
            duration: this.currentEpisode.endTime - this.currentEpisode.startTime,
            numDecisions: this.currentEpisode.decisions.length
        };
        
        this.data.episodes.push(episodeData);
        this.currentRun.episodes.push(episodeData);
        
        // Reset current episode
        this.currentEpisode = {
            decisions: [],
            startTime: null,
            endTime: null,
            startBalance: 0,
            endBalance: 0,
            reward: 0,
            runId: null,
            episodeId: null
        };
        
        this.saveData();
        return episodeData;
    }
    
    finishCurrentRun(endBalance = 0) {
        if (this.currentRun.runId === null) {
            return;
        }
        
        this.currentRun.endTime = Date.now();
        this.currentRun.endBalance = endBalance;
        
        const runData = {
            ...this.currentRun,
            duration: this.currentRun.endTime - this.currentRun.startTime,
            numEpisodes: this.currentRun.episodes.length,
            totalReward: this.currentRun.endBalance - this.currentRun.startBalance
        };
        
        this.data.runs.push(runData);
        
        // Reset current run
        this.currentRun = {
            episodes: [],
            startTime: null,
            endTime: null,
            startBalance: 100,
            endBalance: 0,
            runId: null
        };
        
        this.saveData();
        return runData;
    }
    
    getActionDescription(action) {
        const actionMap = {
            0: 'Respin Blank',
            1: 'Respin Coin',
            2: 'Respin Coin Stack',
            3: 'Respin Snake',
            4: 'Respin Net',
            5: 'Respin x2',
            6: 'Respin Clover',
            7: 'Respin Crown',
            8: 'Cash Out'
        };
        return actionMap[action] || 'Unknown';
    }
    
    getStatistics() {
        const stats = {
            totalDecisions: this.data.decisions.length,
            totalEpisodes: this.data.episodes.length,
            totalRuns: this.data.runs.length,
            avgEpisodeReward: 0,
            avgRunReward: 0,
            actionDistribution: {},
            avgDecisionsPerEpisode: 0,
            avgEpisodesPerRun: 0,
            sessionDuration: Date.now() - this.data.sessionStart
        };
        
        // Calculate average rewards
        if (this.data.episodes.length > 0) {
            const totalEpisodeReward = this.data.episodes.reduce((sum, ep) => sum + ep.reward, 0);
            stats.avgEpisodeReward = totalEpisodeReward / this.data.episodes.length;
            
            const totalDecisions = this.data.episodes.reduce((sum, ep) => sum + ep.numDecisions, 0);
            stats.avgDecisionsPerEpisode = totalDecisions / this.data.episodes.length;
        }
        
        if (this.data.runs.length > 0) {
            const totalRunReward = this.data.runs.reduce((sum, run) => sum + run.totalReward, 0);
            stats.avgRunReward = totalRunReward / this.data.runs.length;
            
            const totalEpisodes = this.data.runs.reduce((sum, run) => sum + run.numEpisodes, 0);
            stats.avgEpisodesPerRun = totalEpisodes / this.data.runs.length;
        }
        
        // Calculate action distribution
        for (const decision of this.data.decisions) {
            const action = decision.action;
            stats.actionDistribution[action] = (stats.actionDistribution[action] || 0) + 1;
        }
        
        return stats;
    }
    
    getPolicyDistribution() {
        const policyMap = new Map();
        
        for (const decision of this.data.decisions) {
            const stateKey = decision.stateCounts.join(',');
            if (!policyMap.has(stateKey)) {
                policyMap.set(stateKey, {
                    count: 0,
                    actions: {}
                });
            }
            
            const stateData = policyMap.get(stateKey);
            stateData.count++;
            stateData.actions[decision.action] = (stateData.actions[decision.action] || 0) + 1;
        }
        
        // Convert to array format for easier processing
        const policyArray = [];
        for (const [stateKey, data] of policyMap.entries()) {
            const stateCounts = stateKey.split(',').map(Number);
            const mostCommonAction = Object.entries(data.actions)
                .reduce((a, b) => data.actions[a[0]] > data.actions[b[0]] ? a : b)[0];
            
            policyArray.push({
                stateCounts: stateCounts,
                stateKey: stateKey,
                totalVisits: data.count,
                actionDistribution: data.actions,
                mostCommonAction: parseInt(mostCommonAction)
            });
        }
        
        return policyArray;
    }
    
    getBalanceHistory() {
        const balanceHistory = [];
        
        for (const episode of this.data.episodes) {
            balanceHistory.push({
                episodeId: episode.episodeId,
                runId: episode.runId,
                balance: episode.endBalance,
                reward: episode.reward,
                timestamp: episode.endTime
            });
        }
        
        return balanceHistory;
    }
    
    exportToJSON() {
        const exportData = {
            metadata: {
                exportTime: Date.now(),
                version: '1.0',
                sessionStart: this.data.sessionStart,
                sessionDuration: Date.now() - this.data.sessionStart
            },
            statistics: this.getStatistics(),
            decisions: this.data.decisions,
            episodes: this.data.episodes,
            runs: this.data.runs,
            policyDistribution: this.getPolicyDistribution(),
            balanceHistory: this.getBalanceHistory()
        };
        
        return JSON.stringify(exportData, null, 2);
    }
    
    exportToCSV() {
        let csv = 'decisionId,runId,episodeId,timestamp,action,actionDescription,respinsUsed,balance,currentPayout,totalCost,';
        csv += 'blank,coin,coinStack,snake,net,x2,clover,crown,';
        csv += 'reel0,reel1,reel2,reel3,validActions\n';
        
        for (const decision of this.data.decisions) {
            const row = [
                decision.decisionId,
                decision.runId,
                decision.episodeId,
                decision.timestamp,
                decision.action,
                decision.actionDescription,
                decision.respinsUsed,
                decision.balance,
                decision.currentPayout,
                decision.totalCost,
                ...decision.stateCounts.slice(0, 8), // Symbol counts
                ...decision.state, // Reel positions
                '"' + decision.validActions.join(';') + '"'
            ];
            csv += row.join(',') + '\n';
        }
        
        return csv;
    }
    
    saveData() {
        try {
            localStorage.setItem('slotMachineData', JSON.stringify(this.data));
            localStorage.setItem('slotMachineCounters', JSON.stringify({
                runCounter: this.runCounter,
                episodeCounter: this.episodeCounter,
                decisionCounter: this.decisionCounter
            }));
        } catch (error) {
            console.error('Error saving data to localStorage:', error);
        }
    }
    
    loadData() {
        try {
            const savedData = localStorage.getItem('slotMachineData');
            if (savedData) {
                const parsedData = JSON.parse(savedData);
                this.data = {
                    decisions: parsedData.decisions || [],
                    episodes: parsedData.episodes || [],
                    runs: parsedData.runs || [],
                    sessionStart: parsedData.sessionStart || Date.now()
                };
            }
            
            const savedCounters = localStorage.getItem('slotMachineCounters');
            if (savedCounters) {
                const parsedCounters = JSON.parse(savedCounters);
                this.runCounter = parsedCounters.runCounter || 0;
                this.episodeCounter = parsedCounters.episodeCounter || 0;
                this.decisionCounter = parsedCounters.decisionCounter || 0;
            }
        } catch (error) {
            console.error('Error loading data from localStorage:', error);
        }
    }
    
    clearData() {
        this.data = {
            decisions: [],
            episodes: [],
            runs: [],
            sessionStart: Date.now()
        };
        this.runCounter = 0;
        this.episodeCounter = 0;
        this.decisionCounter = 0;
        this.currentEpisode = {
            decisions: [],
            startTime: null,
            endTime: null,
            startBalance: 0,
            endBalance: 0,
            reward: 0,
            runId: null,
            episodeId: null
        };
        this.currentRun = {
            episodes: [],
            startTime: null,
            endTime: null,
            startBalance: 100,
            endBalance: 0,
            runId: null
        };
        
        try {
            localStorage.removeItem('slotMachineData');
            localStorage.removeItem('slotMachineCounters');
        } catch (error) {
            console.error('Error clearing data from localStorage:', error);
        }
    }
}