# DSM123 Downstream Modeling and Figures

This directory contains the downstream application and figure-generation layer of Workflow 1. All analyses use the final curated model at `../Models/purple_bacteriav_DSM123.json`.

## Notebooks

- `enzyme_constrained_dfba.ipynb`: enzyme-constrained dFBA, active transport versus diffusion, ForT parameter sensitivity, loopless FVA, and pathway-module capacity scaling.
- `formate_light_scan.ipynb`: steady-state formate-light FBA scan and `RBPC/FDH` flux-ratio heatmap.

Run both notebooks from this directory:

```bash
jupyter nbconvert --to notebook --execute --inplace enzyme_constrained_dfba.ipynb
jupyter nbconvert --to notebook --execute --inplace formate_light_scan.ipynb
```

## Inputs and parameter tables

- `Supplementary File S1. DSM 123 RNA-seq DEG.xlsx`: RNA-seq input used for relative enzyme allocation.
- `DSM123_blast_gene_rename_output_*/`: local locus and protein mapping tables.
- `kinetics.csv`: curated kinetic constants and provenance.
- `substrates.csv`: substrate concentrations used in Michaelis-Menten saturation terms.
- `rnaseq_allocation.csv`: reaction-level RNA-seq allocation values.
- `reaction_parameters.csv`: complete reaction-level kinetic and enzyme-allocation table.
- `dfba_parameters.csv`: compact parameter set used by the active-transport and diffusion comparison.

## Analyses and main figures

| Analysis | Main figure |
|---|---|
| Formate-light FBA scan | `Scan_Analysis_Final.svg` |
| Active ForT transport versus passive diffusion | `mm_enzyme_constrained_dfba_20mM_active_diffusion.svg` |
| ForT `Km` and `Vmax` sensitivity | `mm_dfba_20mM_fort_parameter_sensitivity.svg` |
| Loopless FVA under ForT capacity scaling | `mm_fva_transport_gated_carbon_flux_capacity_loopless_filtered.svg` |
| Module-capacity scaling response | `mm_dfba_module_capacity_scaling_response_heatmap_CS_standard.svg` |
| Module biomass sensitivity coefficients | `mm_dfba_module_biomass_sensitivity_coefficient_heatmap.svg` |

## Module-capacity scaling figures

### Methods

Module-level capacity perturbations were performed with the enzyme-constrained dFBA framework in `enzyme_constrained_dfba.ipynb`. The simulation used the final curated DSM123 model (`../Models/purple_bacteriav_DSM123.json`), an initial formate concentration of 20 mM, active ForT transport with the baseline ForT background fixed at `1.0x`, blocked external CO2 uptake, and a 200 h cultivation window with `dt = 2 h`. For each module, the relevant capacity was scaled over `0.01x`, `0.05x`, `0.1x`, `0.2x`, `0.5x`, and `1.0x` while all other modules were kept at their baseline capacities.

For enzyme modules, capacity scaling was implemented by multiplying the corresponding `kcat` values in the Michaelis-Menten bounds. ForT transport was scaled through the maximal formate uptake capacity. Cyclic electron flow was scaled through the `PURPLE_RC` capacity bound. Photon uptake was scaled through the photon input constraint. The tested modules were `ForT transport`, `FDH`, `CBB`, `TCA`, `Cyclic electron flow`, and `Photon uptake`.

The module-response heatmap reports each output relative to the same module's `1.0x` condition:

```text
relative response = output(alpha) / output(1.0x)
```

The three reported outputs are final biomass, cumulative RuBisCO flux, and cumulative CS flux. The same 200 h result table (`mm_dfba_module_capacity_scaling_data.csv`) was used to calculate the biomass sensitivity coefficient heatmap. Local sensitivity was calculated between adjacent capacity intervals as log-elasticity:

```text
S = [ln(Y_high) - ln(Y_low)] / [ln(alpha_high) - ln(alpha_low)]
```

where `Y` is `final_biomass_gDW_L`. Negative values caused by small numerical plateau fluctuations were clipped to zero before plotting.

### Results

The response heatmap shows that biomass and carbon-assimilation fluxes are most strongly reduced when ForT transport, FDH, CBB, or photon uptake are constrained. FDH remains a strong control point near the upper capacity range, whereas CBB is most influential at intermediate-low capacity. Cyclic electron flow affects the system only under severe capacity knockdown, and TCA scaling does not change the tested outputs under the current model constraints.

The biomass sensitivity coefficient heatmap quantifies these trends. ForT transport shows its strongest local sensitivity in the `0.2-0.5x` interval (`S = 1.66`). FDH is also strongest in the `0.2-0.5x` interval (`S = 1.59`) and remains sensitive in the `0.5-1.0x` interval (`S = 0.95`). CBB is most sensitive across `0.1-0.2x` (`S = 1.65`). Cyclic electron flow is only sensitive at the lowest tested interval (`0.01-0.05x`, `S = 0.34`), while TCA remains insensitive across all tested intervals (`S = 0`).

### Figure Captions

**Module-capacity scaling response under enzyme-constrained dFBA.** Heatmaps show the effect of module-level capacity knockdown on final biomass, cumulative RuBisCO flux, and cumulative citrate synthase (CS) flux during a 200 h enzyme-constrained dFBA simulation. The ForT background was fixed at `1.0x`, external CO2 uptake was blocked, and module capacities were scaled from `0.01x` to `1.0x`. Values are normalized to the corresponding `1.0x` condition for each module. Blue indicates strong loss of output relative to baseline, whereas red indicates retention of baseline output. The analysis identifies ForT transport, FDH, CBB, and photon uptake as the dominant capacity-sensitive modules, while TCA capacity is not limiting under these conditions.

**Local biomass sensitivity coefficient across module-capacity intervals.** Heatmap of log-elasticity coefficients calculated from the same 200 h dFBA module-scaling simulations. Each cell reports the local sensitivity of final biomass to a change in module capacity between adjacent scale intervals using `S = Delta ln(final biomass) / Delta ln(module capacity)`. Larger values indicate stronger biomass control over that capacity interval. The sensitivity analysis shows that ForT transport and FDH exert the strongest control at higher-intermediate capacities, CBB controls biomass primarily at intermediate-low capacity, cyclic electron flow matters only under severe knockdown, and TCA is not limiting in the tested condition.

## Core equations

Active ForT uptake:

```text
v_ForT = Vmax_ForT * S_formate / (Km_ForT + S_formate)
```

Baseline parameters are `Vmax_ForT = 1.786555464 mmol gDW^-1 h^-1` and `Km_ForT = 0.01814 mM`.

RNA-seq-scaled enzyme capacity:

```text
Vmax_r = kcat_r * (enzyme_umol_per_gDW_r / 1000) * 3600
v_bound_r = Vmax_r * min_i[S_i / (Km_i + S_i)]
```

The model biomass protein coefficient is used as the total protein proxy and distributed among constrained reactions using mapped 6 h RNA-seq expression. This is a model-based allocation assumption, not an absolute proteomics measurement.

## Main numerical outputs

- `mm_dfba_20mM_active_diffusion_plot_data.csv`
- `mm_dfba_20mM_fort_parameter_sensitivity_data.csv`
- `mm_fva_transport_gated_capacity_raw_values.csv`
- `mm_fva_transport_gated_capacity_panel_a_values.csv`
- `mm_fva_transport_gated_capacity_normalized_width_values.csv`
- `mm_dfba_module_capacity_scaling_data.csv`

## Interpretation limits

1. `Scan_Analysis_Final.svg` is a conventional steady-state FBA scan; it does not use the RNA-seq-scaled enzyme constraints, and external CO2 uptake is allowed.
2. FVA heatmap rows are normalized independently. A value of 1 is the largest width for that reaction across the tested ForT conditions, not an essentiality score.
3. `CS` is retained as the modeled citrate-synthase entry into the TCA cycle. A broad CS range does not imply that CS is dispensable.
4. Module responses are normalized to the same module's `1x` condition under the same ForT background.
5. The current photon-capacity implementation applies compound scaling proportional to `alpha^2`; interpret that row according to the implemented perturbation.

Detailed methods, equations, parameter provenance, and figure interpretation are provided in `Supplementary_Methods.docx`.

