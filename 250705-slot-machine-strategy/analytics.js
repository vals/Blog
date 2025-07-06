class Analytics {
    constructor(dataCollector) {
        this.dataCollector = dataCollector;
    }
    
    generateBalancePlot() {
        const balanceHistory = this.dataCollector.getBalanceHistory();
        if (balanceHistory.length === 0) {
            return null;
        }
        
        // Create a simple text-based balance plot
        const episodes = balanceHistory.map(entry => entry.episodeId);
        const balances = balanceHistory.map(entry => entry.balance);
        const rewards = balanceHistory.map(entry => entry.reward);
        
        const minBalance = Math.min(...balances);
        const maxBalance = Math.max(...balances);
        const range = maxBalance - minBalance;
        
        let plot = 'BALANCE OVER EPISODES\\n';
        plot += '=' .repeat(50) + '\\n';
        
        for (let i = 0; i < episodes.length; i++) {
            const episode = episodes[i];
            const balance = balances[i];
            const reward = rewards[i];
            
            // Create a simple bar representation
            const normalizedBalance = range > 0 ? (balance - minBalance) / range : 0.5;
            const barLength = Math.round(normalizedBalance * 30);
            const bar = '█'.repeat(barLength) + '░'.repeat(30 - barLength);
            
            plot += `Ep ${episode.toString().padStart(3)}: ${bar} ${balance.toString().padStart(3)} (${reward >= 0 ? '+' : ''}${reward})\\n`;
        }
        
        plot += '\\nBalance Range: ' + minBalance + ' to ' + maxBalance + '\\n';
        plot += 'Total Episodes: ' + episodes.length + '\\n';
        plot += 'Final Balance: ' + balances[balances.length - 1];
        
        return plot;
    }
    
    generatePolicyAnalysis() {
        const policyDistribution = this.dataCollector.getPolicyDistribution();
        const stats = this.dataCollector.getStatistics();
        
        if (policyDistribution.length === 0) {
            return 'No policy data available yet.';
        }
        
        let analysis = 'POLICY ANALYSIS\\n';
        analysis += '=' .repeat(50) + '\\n';
        
        // Overall action distribution
        analysis += 'ACTION DISTRIBUTION:\\n';
        for (const [action, count] of Object.entries(stats.actionDistribution)) {
            const percentage = (count / stats.totalDecisions * 100).toFixed(1);
            const actionName = this.dataCollector.getActionDescription(action);
            analysis += `  ${actionName}: ${count} (${percentage}%)\\n`;
        }
        
        analysis += '\\nSTATE-ACTION PATTERNS:\\n';
        
        // Sort states by visit frequency
        const sortedStates = policyDistribution.sort((a, b) => b.totalVisits - a.totalVisits);
        
        // Show top 10 most visited states
        const topStates = sortedStates.slice(0, 10);
        analysis += 'Top 10 Most Visited States:\\n';
        
        for (const state of topStates) {
            const symbolNames = ['Blank', 'Coin', 'CoinStack', 'Snake', 'Net', 'x2', 'Clover', 'Crown'];
            const stateDesc = state.stateCounts.slice(0, 8).map((count, i) => 
                count > 0 ? `${symbolNames[i]}:${count}` : '').filter(s => s).join(', ');
            const respins = state.stateCounts[8];
            
            analysis += `  [${stateDesc}] (${respins} respins) - ${state.totalVisits} visits\\n`;
            analysis += `    Primary action: ${this.dataCollector.getActionDescription(state.mostCommonAction)}\\n`;
            
            // Show action distribution for this state
            const actionDist = Object.entries(state.actionDistribution)
                .map(([action, count]) => `${this.dataCollector.getActionDescription(action)}:${count}`)
                .join(', ');
            analysis += `    Actions: ${actionDist}\\n`;
        }
        
        // Analysis of decision patterns
        analysis += '\\nDECISION PATTERNS:\\n';
        
        const cashOutStates = policyDistribution.filter(s => s.mostCommonAction === 8);
        const respinStates = policyDistribution.filter(s => s.mostCommonAction !== 8);
        
        analysis += `States where cash-out is preferred: ${cashOutStates.length}\\n`;
        analysis += `States where respins are preferred: ${respinStates.length}\\n`;
        
        // Find patterns in cash-out decisions
        if (cashOutStates.length > 0) {
            const avgCashOutPayout = this.calculateAveragePayoutForStates(cashOutStates);
            analysis += `Average payout in cash-out states: ${avgCashOutPayout.toFixed(2)}\\n`;
        }
        
        return analysis;
    }
    
    calculateAveragePayoutForStates(states) {
        // This would need actual payout data to be meaningful
        // For now, return a placeholder
        return 0;
    }
    
    generateSessionSummary() {
        const stats = this.dataCollector.getStatistics();
        const sessionDurationMinutes = Math.round(stats.sessionDuration / 60000);
        
        let summary = 'SESSION SUMMARY\\n';
        summary += '=' .repeat(50) + '\\n';
        summary += `Session Duration: ${sessionDurationMinutes} minutes\\n`;
        summary += `Total Decisions: ${stats.totalDecisions}\\n`;
        summary += `Total Episodes: ${stats.totalEpisodes}\\n`;
        summary += `Total Runs: ${stats.totalRuns}\\n`;
        summary += `Average Episode Reward: ${stats.avgEpisodeReward.toFixed(2)}\\n`;
        summary += `Average Run Reward: ${stats.avgRunReward.toFixed(2)}\\n`;
        summary += `Average Decisions per Episode: ${stats.avgDecisionsPerEpisode.toFixed(2)}\\n`;
        summary += `Average Episodes per Run: ${stats.avgEpisodesPerRun.toFixed(2)}\\n`;
        
        if (stats.totalDecisions > 0) {
            const decisionsPerMinute = stats.totalDecisions / sessionDurationMinutes;
            summary += `Decision Rate: ${decisionsPerMinute.toFixed(2)} decisions/minute\\n`;
        }
        
        return summary;
    }
    
    generateRecommendations() {
        const stats = this.dataCollector.getStatistics();
        const policyDistribution = this.dataCollector.getPolicyDistribution();
        
        let recommendations = 'RECOMMENDATIONS & INSIGHTS\\n';
        recommendations += '=' .repeat(50) + '\\n';
        
        if (stats.totalDecisions < 50) {
            recommendations += 'Play more episodes to gather meaningful data for analysis.\\n';
            return recommendations;
        }
        
        // Analyze action preferences
        const actionCounts = stats.actionDistribution;
        const totalActions = stats.totalDecisions;
        
        const cashOutPercentage = (actionCounts[8] || 0) / totalActions * 100;
        const respinPercentage = 100 - cashOutPercentage;
        
        recommendations += `Cash-out rate: ${cashOutPercentage.toFixed(1)}%\\n`;
        recommendations += `Respin rate: ${respinPercentage.toFixed(1)}%\\n`;
        
        if (cashOutPercentage > 80) {
            recommendations += 'You tend to cash out frequently. Consider taking more calculated risks.\\n';
        } else if (cashOutPercentage < 20) {
            recommendations += 'You take many risks with respins. Consider cashing out more often.\\n';
        }
        
        // Find most and least used respin actions
        const respinActions = Object.entries(actionCounts)
            .filter(([action, count]) => parseInt(action) < 8)
            .sort((a, b) => b[1] - a[1]);
        
        if (respinActions.length > 0) {
            const mostUsedRespin = respinActions[0];
            const leastUsedRespin = respinActions[respinActions.length - 1];
            
            recommendations += `Most used respin: ${this.dataCollector.getActionDescription(mostUsedRespin[0])} (${mostUsedRespin[1]} times)\\n`;
            recommendations += `Least used respin: ${this.dataCollector.getActionDescription(leastUsedRespin[0])} (${leastUsedRespin[1]} times)\\n`;
        }
        
        // Analyze episode rewards
        if (stats.avgEpisodeReward > 0) {
            recommendations += `Good strategy! Average episode reward is positive: ${stats.avgEpisodeReward.toFixed(2)}\\n`;
        } else {
            recommendations += `Strategy needs improvement. Average episode reward is negative: ${stats.avgEpisodeReward.toFixed(2)}\\n`;
        }
        
        return recommendations;
    }
    
    exportAnalysisReport() {
        const report = {
            timestamp: new Date().toISOString(),
            sessionSummary: this.generateSessionSummary(),
            balancePlot: this.generateBalancePlot(),
            policyAnalysis: this.generatePolicyAnalysis(),
            recommendations: this.generateRecommendations(),
            rawStats: this.dataCollector.getStatistics(),
            policyDistribution: this.dataCollector.getPolicyDistribution(),
            balanceHistory: this.dataCollector.getBalanceHistory()
        };
        
        return JSON.stringify(report, null, 2);
    }
    
    displayAnalysisInAlert() {
        const summary = this.generateSessionSummary();
        const recommendations = this.generateRecommendations();
        
        const message = summary + '\\n\\n' + recommendations;
        alert(message);
    }
    
    // Helper method to create downloadable analysis report
    downloadAnalysisReport() {
        const report = this.exportAnalysisReport();
        const blob = new Blob([report], { type: 'application/json' });
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = 'slot_machine_analysis.json';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
    }
}