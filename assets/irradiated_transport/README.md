# Irradiated transport data

- `decay_gamma_reference.csv`: normalized 48-group cooldown photon-spectrum
  shapes and lead mass-attenuation coefficients used by the
  irradiated-transport photon model.
- `shield_material_properties.csv`: shielding density and raw material cost. Lead is the only enabled material in the first implementation. The 2024 reference cost is escalated to 2025.
- `irradiated_transport_cost_inputs_2025.csv`: radioactive-content-specific road and sea screening allowances, with rail values retained for traceability but not separately added to the broad existing rail range.

No file under `assets/Ref_Results/BOC/` is required by the irradiated-transport
runtime calculation.

See `IRRADIATED_TRANSPORT_MODE.md` at the repository root for assumptions, calculation flow, and exclusions.
