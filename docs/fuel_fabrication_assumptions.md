# Fuel fabrication assumptions

Account 254 scales fabrication cost by **kilograms of uranium**, not total fuel,
TRISO-particle, compact, or assembly mass. Uranium purchase, conversion, and
enrichment remain separate cost accounts. The workbook
`cost/Cost_Database.xlsx` is authoritative; its CSV files are generated mirrors.

## Adopted inputs

| Fuel | Dollar year | Low ($/kgU) | Mode ($/kgU) | High ($/kgU) | Distribution | Learning per doubling |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| UO2 | 2017 | 230 | 400 | 575 | Triangular | 2% |
| TRISO | 2017 | 3,300 | 10,900 | 29,400 | Triangular | 13% |

For a deterministic run (`Number of Samples = 1`), the engine uses the mode.
For uncertainty runs, the low/mode/high values parameterize a bounded triangular
distribution. They are not lognormal quantiles. The expected values before
escalation are approximately $401.67/kgU for UO2 and $14,533.33/kgU for TRISO;
the mode is not the arithmetic mean. Existing scaling-exponent uncertainty is
unchanged. TRIGA fabrication inputs and its learning category are unchanged.

The dollar year must travel with the price. MOUSE escalates these inputs using
the workbook's **Material** index. At the default 2025 escalation year, the
deterministic modes are approximately $610.70/kgU (UO2) and $16,641.61/kgU
(TRISO), before deployment learning.

## Published price and uncertainty basis

- **UO2:** [AFC-CBR Module D1-1](https://inldigitallibrary.inl.gov/content/uploads/50/2026/04/Sort_90170.pdf),
  Table D1-1.21, standard PWR UOX: triangular low/mode/high of 230/400/575 in
  2017 USD/kgU. This commercial LWR-fuel benchmark is used as a scoping proxy;
  it does not establish a microreactor-specific or HALEU fabrication quote.
- **TRISO:** [AFC-CBR Module D1-3, September 2021](https://sai.inl.gov/content/uploads/29/2024/11/2021_Module_D1-3.pdf),
  printed page D1-3-27, Table D1-3-10: historical HTR low/mode/high of
  3,300/10,900/29,400 in 2017 USD/kgU. This replaces the shielding branch's
  5,450/10,900/32,700 inputs with the published historical triangular inputs.
  The report also presents a separate 2021 update in **2020 dollars** on page
  D1-3-28. The adopted values deliberately retain the historical 2017 benchmark;
  they are not described as the report's newest estimate or a current vendor
  price. Production scale and whether final compacts, blocks, or pebbles are
  included must be matched before substituting another source.

## Deployment learning is an additional scenario assumption

Both rates below are **inferred effective assumptions**, not measured
cost-versus-cumulative-reactor-deployment relationships published by AFC-CBR.

- **UO2, 2%:** inspired by the throughput-related reductions discussed in
  Module D1-1, Section D1-1.10.2 and Figure D1-1.6. At unit 100, this produces
  a multiplier of about 0.874 (12.6% below the unit-1 cost). Plant throughput
  economies are being approximated by deployment learning, not equated to a
  directly observed learning curve.
- **TRISO, 13%:** inspired by the historical current-production versus automated
  production discussion in Module D1-3, Section D1-3-8. At unit 100, this
  produces a multiplier of about 0.396 (60.4% below the unit-1 cost). The
  historical comparison of roughly $25,000/kgU with a $10,000/kgU automation
  target motivates the fractional reduction, but does not independently
  justify applying it again to a $10,900/kgU starting benchmark.

In particular, the historical TRISO benchmark already reflects expectations
of mature or automated production. Applying another 13% per doubling can
double-count that improvement. The adopted default explicitly models **further**
improvement beyond that benchmark; it should be treated as a sensitivity
assumption, not an established FOAK-to-NOAK forecast. Use `TRISO Learning = 0`
as a no-additional-learning comparison. UO2 learning can likewise be set to zero.

The engine uses `multiplier = (1 - rate) ** log2(min(100, unit_number))`.
Unit 1 has multiplier 1, and improvements are capped at unit 100. The rates
live in the workbook's Economics Parameters sheet, with matching web-app
defaults. Explicit caller overrides take precedence. They apply only to the
fabrication rows assigned to these categories, not the entire fuel cost or
the existing generic Factory Other category.

At unit 100 and in 2025 dollars, the adopted deterministic fabrication modes
are approximately $533.99/kgU (UO2) and $6,597.36/kgU (TRISO). These figures are
account-level rates, not total core costs or total LCOE.
