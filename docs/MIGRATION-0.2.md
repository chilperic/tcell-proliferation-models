# Scientific and API migration to 0.2.0

Cyton defaults to an ODE representation of the exponential competing-risk model, with separate proliferating and nondividing compartments. This removes cancellation and equal-rate singularities in the historical analytic expression. method="legacy" retains that expression for comparison. Pure-division and no-division limits are tested. The new formulation may change historical fitted results; refit and compare explicitly.

Smith-Martin retains constant pre-zero history by default for compatibility. history="zero" represents no pre-zero A-phase cells and prevents delayed division before delta for synchronized initial conditions. Always record the history and legacy_survival choice. No history convention is universally appropriate.

Input time order is checked rather than sorted; replicate counts require identical time and generation coordinates. Missing observations and weight recycling are rejected. Fits are selected only from converged finite optimization results. Full thesis figure reproduction and empirical refits remain pending.
