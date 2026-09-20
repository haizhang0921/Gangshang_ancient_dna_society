Gangshang Northern Cemetery — three kinship-scenario chronology framework
Date: 2026-09-20

Purpose
The same updated OxCal archaeological/radiocarbon input (gangshang_n.csv; 48 northern human 14C determinations, 34 top-level burial events) is combined with three alternative genetic relationship specifications.

Scenario 1: main (genetic_relations_n)
Used in the main manuscript.
20 edges. Explicitly directed: NM30->NM27, NM14->NM32, NM15->NM16.
NM22-NM23, NM31-NM14, NM34-NM14 are other 2nd-degree relationships (35 +/- 26 y).

Scenario 2: alt_close (genetic_relations_n_alt)
Sensitivity model representing genetically closer relationships and therefore potentially smaller death-date separations.
Differences from main:
- NM22-NM23: other 2nd-degree (35 +/- 26) -> sibling (26 +/- 22)
- NM31-NM14: other 2nd-degree (35 +/- 26) -> sibling (26 +/- 22)
- NM34-NM14: other 2nd-degree (35 +/- 26) -> sibling (26 +/- 22)
All directions unchanged.

Scenario 3: alt_far (genetic_relations_n_alt2)
Sensitivity model representing a more distant / less chronologically restrictive interpretation.
Differences from main:
- NM22-NM23: other 2nd-degree (35 +/- 26) -> grandparent-grandchild (35 +/- 32)
- NM30-NM27: explicit death-order direction removed; parent-offspring gap remains 29 +/- 19.
Other edges unchanged.

Recommended reporting
- Main manuscript: Scenario 1.
- Supplementary sensitivity analysis: Scenarios 2 and 3, emphasizing that they bracket plausible kinship uncertainty.
- Use identical OxCal input, likelihood parameter source, MCMC settings, stage-order sensitivity framework and boundary reconstruction method across all scenarios.

Formal R workflow
1. source("run_three_scenarios.R")
2. source("04_compare_three_scenarios.R")
Formal defaults remain 4 chains x 80,000 iterations, burn-in 20,000, thin 20, seed 20260806.

Reference outputs
The reference_python folder contains mathematically equivalent discrete-Gibbs calculations for checking the direction and magnitude of scenario effects. These are not a substitute for the formal R outputs used in the manuscript.
