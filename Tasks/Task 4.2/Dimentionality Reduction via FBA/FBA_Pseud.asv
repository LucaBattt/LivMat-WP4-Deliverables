% Load the Pseudomonas model
model_ps = readCbModel('P_taiwanensis.xlsx');

% Lets list the exchange reactions (exchanges that happen between the
% internal metabolic network and the environment
fprintf("------------------------------------------------------------------\n")
fprintf("Exchange reactions:")
ex = findExcRxns(model_ps);
model_ps.rxns(ex) % A total of 379 exchange reactions

% Make new instance of the model to leave the original untouched
model_psA = model_ps;

% Change the objective to Adipic acid
model_psA = changeObjective(model_psA, 'EX_Adipate');

% Set oxygen uptake
model_psA = changeRxnBounds(model_psA, 'EX_o2(e)', -10, 'l');

% I will make a loop over glucose uptake bounds and save the EX_Adipate
% value

% Glucose uptake values
glucoseRange = -linspace(0,20,41); % 41 points from 0 to -20

% For saving the results
biomassFlux = zeros(length(glucoseRange),1);
adipateFlux = zeros(length(glucoseRange),1);

% Run the loop
for i = 1:length(glucoseRange)
    glc_uptake = glucoseRange(i);

    % Set glucose lower bound
    modelTmp = changeRxnBounds(model_psA, 'EX_glc(e)', glc_uptake, 'l');

    % Solve FBA
    sol = optimizeCbModel(modelTmp);

    % Save the optimized adipate and biomass flux values
    adipateFlux(i) = sol.f; % Adipate flux
    biomassFlux(i) = sol.x(findRxnIDs(modelTmp, 'BiomassKT2440_WT3')); % Biomass flux

end

% Plot adipate and biomass vs glucoseRange
figure;
plot(glucoseRange, adipateFlux, 'r-', 'LineWidth', 2);
hold on;
plot(glucoseRange, biomassFlux, 'b-', 'LineWidth', 2);
xlabel('Glucose Uptake (mmol/gDW/h)');
ylabel('Flux (mmol/gDW/h)');
legend('Adipate Flux', 'Biomass Flux');
title('Adipate and Biomass Flux vs Glucose Uptake');
grid on;


% --------------- Performing FBA to select top n fluxes -------------

% Choose a single value for glucose uptake
glc_uptake = -10;

% Run FBA
model_ps10 = changeRxnBounds(model_psA, 'EX_glc(e)', glc_uptake, 'l');
sol_ps10   = optimizeCbModel(model_ps10);
fluxes     = sol_ps10.x;

% Extract top n fluxes by absolute value
n = 100;
[~, idx] = sort(abs(fluxes), 'descend'); % Sort fluxes by absolute value
topFluxes = fluxes(idx(1:n)); % Select top 20 fluxes
topRxns = model_ps10.rxns(idx(1:n)); % Get corresponding reaction IDs

% Plot histogram
figure;
bar(abs(topFluxes));
set(gca,'XTickLabel',topRxns,'XTickLabelRotation',90);
ylabel('flux value');
title('Top 100 fluxes (abs)');


% --------------------- Perform reaction KO -------------------------
% Knock-out each of the top 100 reactions and re-optimize.
% Record how much the objective drops.
% Establish a hierarchical order of importance based on that.

impact = zeros(n,1);

for i = 1:n
    rxn = topRxns{i};

    % Perform KO
    model_psKO = changeRxnBounds(model_ps10, rxn, 0, 'b'); % Knockout reaction
    sol_KO = optimizeCbModel(model_psKO); % Re-optimize the model
    
    if isempty(sol_KO.f)
        impact(i) = sol_ps10.f; % Infeasible
    else
        impact(i) = sol_ps10.f - sol_KO.f; % Record the change in objective value
    end  
end


% Display the impact of each knockout reaction
figure;
bar(impact);
set(gca,'XTickLabel',topRxns,'XTickLabelRotation',90);
xlabel('KO Reaction');
ylabel('Impact on Objective Value');
title('Impact of Reaction Knockouts on Objective Value');
grid on;

% Keep only the fluxes and reactions that make the objective value drop by
% more than 10% of its value when it is knocked out.

threshold = 0.1 * sol_ps10.f;
significantKnockouts = topRxns(impact > threshold);
significantImpacts = impact(impact > threshold);

% Display significant knockouts and their impacts
disp('Significant Knockouts:');
disp(significantKnockouts);
disp('Impacts on Objective Value:');
disp(significantImpacts);

% Display sorted results in a table
fprintf("Table of top 100 reactions\n")
T = table(topRxns(:), impact(:), 'VariableNames', {'Reaction', 'Flux'});
[~, idx] = sort(impact, 'descend');
Tsorted = T(idx,:);
disp(Tsorted)
writetable(Tsorted, 'top_reactions_ps.csv')

% Display sorted significant KO results in a table
fprintf("Table of significant KO results\n")
T_KO = table(significantKnockouts(:), significantImpacts(:), 'VariableNames', {'Reaction', 'Impact'});
[~, idx] = sort(significantImpacts, 'descend');
T_KO_sorted = T_KO(idx,:);
disp(T_KO_sorted)
writetable(T_KO_sorted, 'top_reactions_ps_ko.csv')

% Print total number of reactions in the model
fprintf("Reactions:             %d\n", length(model_ps.rxns));
fprintf("Genes:                 %d\n", length(model_ps.genes));

% Print total number of identified significant reactions
fprintf("Significant Reactions: %d\n", length(significantKnockouts))

% Save the significant knockouts and their impacts to a CSV file
%writetable(Tsorted, 'significant_knockouts.csv');