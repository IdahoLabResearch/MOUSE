# Depletion-derived fuel lifetime

MOUSE reports the method used to estimate operating lifetime from its depletion
history. The same lifetime rules apply to the corrected 2D and explicit 3D
operating-depletion workflows.

| Depletion outcome | Reported lifetime | Cost estimation |
| --- | --- | --- |
| Starting operating keff is at or below 1 | Zero days; `BOL subcritical` | Blocked with a clear error |
| Simulated operating keff crosses downward through 1 | Linear interpolation between the bracketing points; `Interpolated k=1 crossing` | Allowed for a positive lifetime |
| No crossing, but the final two points decrease | Linear prediction from the final slope; `Extrapolated beyond depletion schedule` | Blocked; extend depletion until a crossing is simulated |
| No crossing and no decreasing terminal slope | Last simulated time as a lower bound; `Simulated endpoint lower bound` | Blocked; extend depletion until a crossing is simulated |

`Fuel Lifetime` remains in days. Corrected 2D results retain their existing
nearest-day precision; explicit 3D results retain fractional days. There is no
rounding into calendar-year study bins. The additional
fields are `Fuel Lifetime Status`, `Fuel Lifetime Extrapolated`, and
`Fuel Lifetime Simulated Endpoint (days)`. The endpoint is the last simulated
time and is retained independently of an extrapolated prediction. For example,
an operating history that remains above 1 through day 2,000 may predict a
crossing at day 2,500, but day 2,000 remains the simulated endpoint.

Extrapolation is a provisional report based only on the final two values. It
does not simulate the future depleted fuel composition, propagate statistical
uncertainty into a lifetime uncertainty, or establish end-of-life conditions.
A shallow terminal slope may produce a distant prediction; that prediction
still cannot be used for costs.

The operating workflow records the result and displays the reason before
raising an error for a subcritical or unbracketed result. It stops before
downstream lifetime-dependent calculations. In particular, shutdown-margin
and temperature-coefficient lifecycle evaluations require a simulated
end-of-life crossing; extrapolation does not supply the missing composition.
The error asks the user to extend the depletion schedule for unbracketed
results. A subcritical start instead requires revising the operating design.

The shared calculation retains its unrounded estimate internally, and the
simulated endpoint is stored in days. The explicit 3D shutdown-state pass is
allowed to be subcritical: it supplies shutdown-margin data, not an operating
lifetime for costs. Lifetime reporting and rejection apply to the operating
pass, which follows it.

The reactor and enabled central-facility cost entry points independently
reject zero, negative, nonfinite, or explicitly provisional lifetime inputs
before reading the database, performing cost calculations, or writing cost
outputs. This also protects callers that catch a workflow error and then try
to run costs. A disabled central-facility calculation remains a no-op.

Existing positive lifetime inputs without depletion-status metadata remain
supported, including the web app's reference-data estimation path. Its
existing subcritical/short-lifetime checks remain in effect. This change does
not label those reference-data estimates as simulated depletion crossings.

No study-constraint summary flags, rounded-cycle acceptance ranges, heat-pipe
loading limits, or depletion-integrator changes are included.
