% Load the Synechocystis model
model_sy = readCbModel('Synechocystis.xlsx');

% Lets list the exchange reactions (exchanges that happen between the
% internal metabolic network and the environment
fprintf("------------------------------------------------------------------\n")
fprintf("Exchange reactions:")
ex = findExcRxns(model_sy);
model_sy.rxns(ex) % A total of 58 exchange reactions

% Make new instance of the model to leave the original untouched
model_syA = model_sy;

% Change the objective to biomass
model_syA = changeObjective(model_syA, 'Ec_biomass_SynMixo');

% Set oxygen uptake
model_syA = changeRxnBounds(model_syA, 'EX_o2(e)', -10, 'l');

% Set glucos uptake
model_syA = changeRxnBounds(model_syA, 'EX_glc(e)', -10, 'l');

% I will make a loop over glucose uptake bounds and save the biomass value

% Glucose uptake values
glucoseRange = -linspace(0,20,41); % 41 points from 0 to -20

% For saving the results
biomassFlux = zeros(length(glucoseRange),1);
adipateFlux = zeros(length(glucoseRange),1);

% Run the loop
for i = 1:length(glucoseRange)
    glc_uptake = glucoseRange(i);

    % Set glucose lower bound
    modelTmp = changeRxnBounds(model_syA, 'EX_glc(e)', glc_uptake, 'l');

    % Solve FBA
    sol = optimizeCbModel(modelTmp);

    % Save the optimized adipate and biomass flux values
    biomassFlux(i) = sol.f; % Biomass flux
end

% Plot biomass vs glucoseRange
figure;
plot(glucoseRange, biomassFlux, 'b-', 'LineWidth', 2);
xlabel('Glucose Uptake (mmol/gDW/h)');
ylabel('Flux (mmol/gDW/h)');
legend('Biomass Flux');
title('Biomass Flux vs Glucose Uptake');
grid on;


% --------------- Performing FBA to select top n fluxes -------------

% Choose a single value for glucose uptake
glc_uptake = -10;

% Run FBA
model_sy10 = changeRxnBounds(model_syA, 'EX_glc(e)', glc_uptake, 'l');
sol_sy10 = optimizeCbModel(model_sy10);
fluxes     = sol_sy10.x;


% Extract top n fluxes by absolute value
n = 100;
[~, idx] = sort(abs(fluxes), 'descend'); % Sort fluxes by absolute value
topFluxes = fluxes(idx(1:n)); % Select top n fluxes
topRxns = model_sy10.rxns(idx(1:n)); % Get corresponding reaction IDs

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
    model_syKO = changeRxnBounds(model_sy10, rxn, 0, 'b'); % Knockout reaction
    sol_KO = optimizeCbModel(model_syKO); % Re-optimize the model
    
    if isempty(sol_KO.f)
        impact(i) = sol_sy10.f; % Infeasible
    else
        impact(i) = sol_sy10.f - sol_KO.f; % Record the change in objective value
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

threshold = 0.1 * sol_sy10.f;
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
writetable(Tsorted, 'top_reactions_sy.csv')

% Display sorted significant KO results in a table
fprintf("Table of significant KO results\n")
T_KO = table(significantKnockouts(:), significantImpacts(:), 'VariableNames', {'Reaction', 'Impact'});
[~, idx] = sort(significantImpacts, 'descend');
T_KO_sorted = T_KO(idx,:);
disp(T_KO_sorted)
writetable(T_KO_sorted, 'top_reactions_sy_ko.csv')

% Print total number of reactions in the model
fprintf("Reactions:             %d\n", length(model_sy.rxns));
fprintf("Genes:                 %d\n", length(model_sy.genes));

% Print total number of identified significant reactions
fprintf("Significant Reactions: %d\n", length(significantKnockouts))

% Save the significant knockouts and their impacts to a CSV file
%writetable(Tsorted, 'significant_knockouts.csv');