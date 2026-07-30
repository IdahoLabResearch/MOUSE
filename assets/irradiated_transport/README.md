# Irradiated transport data

- `decay_gamma_reference.csv`: normalized 48-group cooldown photon-spectrum shapes and lead mass-attenuation coefficients derived from the Manit Shah workbook.
- `shield_material_properties.csv`: shielding density and raw material cost. Lead is the only enabled material in the first implementation. The 2024 Shah cost is escalated to 2025.
- `irradiated_transport_cost_inputs_2025.csv`: radioactive-content-specific road and sea screening allowances, with rail values retained for traceability but not separately added to the broad existing rail range.

See `IRRADIATED_TRANSPORT_MODE.md` at the repository root for assumptions, calculation flow, and exclusions.
