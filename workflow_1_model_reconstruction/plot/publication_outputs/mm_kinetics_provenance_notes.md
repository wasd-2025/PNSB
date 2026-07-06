# Michaelis-Menten enzyme constraint notes

Each listed CBB/TCA reaction is constrained at every dFBA step with `v = kcat_s_inv * (enzyme_umol_per_gdw / 1000) * 3600 * min(S_i/(Km_i+S_i))`. The model has ACXYSJ locus tags and local GenBank-style protein IDs. Exact UniProt accessions were not present in the local annotation, so the exported table keeps a UniProt query field for each reaction. Replace kcat, Km, enzyme amount, and intracellular substrate assumptions with accession-specific UniProt/BRENDA/SABIO-RK values when those are confirmed.
