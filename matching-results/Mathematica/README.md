# MSSM-to-SMEFT matching results in Mathematica format
This directory contains all matching results in the format of Wolfram language files with Matchete syntax. The following files are included:

- `MSSM-matching-conditions.m`: full one-loop matching conditions as derived in the notebook [MSSM-matching.nb](../../matching/MSSM-matching.nb), where the renormalizable SMEFT couplings have been shifted to absorb their BSM contributions, and comparison to the file listed next;
- `MSSM-matching-conditions_detailed.m`: full one-loop matching conditions as derived step-by-step in the notebook [MSSM-matching-detailed.nb](../../matching/MSSM-matching-detailed.nb), where again the renormalizable SMEFT couplings have been shifted to absorb their BSM contributions;
- `MSSM-matching-conditions_detailed_not-shifted.m`: full one-loop matching conditions as derived step-by-step in the notebook [MSSM-matching-detailed.nb](../../matching/MSSM-matching-detailed.nb), but without shifting the renormalizable SMEFT couplings;
- `SMEFT-Lagrangian_off-shell.m`: full SMEFT Lagrangian obtained from the matching before *any* simplifications are applied, as derived in the notebook [MSSM-matching.nb](../../matching/MSSM-matching.nb);
- `SMEFT-Lagrangian_on-shell.m`: full SMEFT Lagrangian obtained from the matching reduced to the minimal on-shell basis employed by Matchete.