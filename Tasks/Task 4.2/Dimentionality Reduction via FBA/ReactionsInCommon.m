% Load models
model_ps = readCbModel('P_taiwanensis.xlsx');
model_sy = readCbModel('Synechocystis.xlsx');

% Extract reaction IDs
rx_ps = model_ps.rxns;
rx_sy = model_sy.rxns;

% Common reactions
common_rxns = intersect(rx_ps, rx_sy);

% Unique to each model
only_ps = setdiff(rx_ps, rx_sy);
only_sy = setdiff(rx_sy, rx_ps);

% Display results
fprintf('Number of reactions in P. taiwanensis: %d\n', numel(rx_ps));
fprintf('Number of reactions in Synechocystis: %d\n', numel(rx_sy));
fprintf('Number of common reactions: %d\n', numel(common_rxns));
fprintf('Reactions only in P. taiwanensis: %d\n', numel(only_ps));
fprintf('Reactions only in Synechocystis: %d\n', numel(only_sy));
fprintf('Total unique reactions: %d\n', numel(only_ps) + numel(only_sy));