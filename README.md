<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./assets/logos/MOUSE-logo_R1_white.png">
    <source media="(prefers-color-scheme: light)" srcset="./assets/logos/MOUSE-logo_R1_black.png">
    <img alt="MOUSE logo" src="./assets/logos/MOUSE-logo_R1_black.png" width="500">
  </picture>
</p>

# Microreactor Optimization Using Simulation and Economics (MOUSE) 

## Web App

Use the hosted MOUSE Streamlit app here:

**[Launch the MOUSE Streamlit App](https://mouse-microreactors.streamlit.app)**

Project overview page:
**[https://idaholabresearch.github.io/MOUSE/](https://idaholabresearch.github.io/MOUSE/)**

The web app is an interactive scoping tool that bridges microreactor design and
economics without requiring users to install OpenMC or WATTS locally. It connects
early-stage design choices with bottom-up first-of-a-kind (FOAK) and nth-of-a-kind
(NOAK) economics for rapid technical and economic screening.

## Motivation

Nuclear microreactors are being developed for transportable, rapidly deployable,
resilient power in remote or infrastructure-limited locations. However, their
small size and modularity can exacerbate diseconomies of scale, and
first-of-a-kind costs are often driven by expensive components, specialty
materials, fuel fabrication, and indirect costs that do not shrink
proportionally with power.

Historically, reactor designs have often been guided by physics and engineering
first, with economics assessed later using simplified scaling. For
microreactors, that approach is risky because design choices can significantly
affect component mass, fuel lifetime, balance-of-plant needs, capital cost, and
levelized cost of energy (LCOE). MOUSE addresses this need by tightly coupling
design calculations with bottom-up cost estimation.

## Description

MOUSE integrates nuclear microreactor design and economics. It leverages the
[OpenMC](https://github.com/openmc-dev/openmc) Monte Carlo Particle Transport
Code for neutronics core simulations and uses the Workflow and Template Toolkit
for Simulation ([WATTS](https://github.com/watts-dev/watts)) for parametric
studies. MOUSE includes simplified balance-of-plant and operational
calculations, and its cost framework draws on data from the MARVEL project and
other open literature.

Economically, MOUSE provides comprehensive bottom-up cost estimates for
near-term and longer-term microreactors, including capital costs, annualized
costs, fuel costs, and LCOE. It supports analysis of technological factors,
design changes, materials, geometry modifications, and economic parameters.

## What MOUSE Does

MOUSE supports:
- Parametric microreactor design studies using LTMR, GCMR, and HPMR reference concepts
- Neutronics-informed core simulations and first-order thermal-hydraulic scoping
- Balance-of-plant calculations and supporting plant equipment estimates
- Bottom-up capital, annualized, fuel, and operating cost estimation
- FOAK and NOAK calculations for OCC, TCI, LCOE, LCOH, and LCOF
- Cost-driver analysis, uncertainty ranges, IRA tax credit effects, and market comparisons
- Transportability screening for component dimensions, mass, truck, rail, and sea movement
- LTMR, GCMR, and HPMR reference designs
- Interactive web-based scoping through the Streamlit app

## MOUSE Tool Structure
<img src="./assets/mouse_diagram.png" />

### User-Defined Inputs
Users can modify design inputs or economic inputs such as:
- Overall System: Reactor Power (MWt), Thermal Efficiency (%), Heat Flux Criteria
- Geometry: Fuel Pin Radii, TRISO Packing Fraction, Coolant Channel Radius, Moderator Booster Radius, Lattice Pitch, Rings per Assembly, Assemblies per Core, Core Active Height, Reflector Thickness, Control Drum Dimensions
- Materials: Fuel, Enrichment, Coolant, Reflector, Matrix Material, Moderator, Moderator Booster, Control Drum Absorber/Reflector, Fuel Pin Materials
- Shielding: In/Out Vessel Shield Thickness, Material, Dimensions
- Vessels: Vessel Radius, Thickness, Materials, Gaps Between Vessels
- Balance of Plant: Coolant Inlet/Outlet Temperatures, Compressor Pressure Ratio, Pump Efficiency
- Operation: Operation Mode, Number of On-Site Operators per Shift, Plant Lifetime, Refueling Period, Number of Emergency Shutdowns per Year, Startup Durations
- Buildings: Dimensions of Reactor, Turbine, Control, Refueling, Spent Fuel, Emergency, Storage, Radioactive Waste Buildings
- Economic Parameters: Interest Rate, Dollar Escalation Year, Construction Duration, Debt to Equity Ratio

MOUSE is powered by the Workflow and Template Toolkit for Simulation ([WATTS](https://github.com/watts-dev/watts)), developed by ANL, which facilitates parametric studies by integrating various code components.

### Microreactor Designs
Three reactor designs are included so far:
- Liquid-metal thermal microreactor (LTMR)
- Gas-cooled microreactor (GCMR)
- Heat pipe microreactor (HPMR)

The designs can be found [here](./assets/Ref_openmc_2d_designs)


## Ways to Use MOUSE

### Option 1: Use the Web App

Use the hosted Streamlit app for interactive microreactor scoping and cost estimation:

**[https://mouse-microreactors.streamlit.app](https://mouse-microreactors.streamlit.app)**

No local installation is required.

### Option 2: Run the Research Code Locally

For full OpenMC/WATTS-based workflows, install the required dependencies listed below.

## Prerequisites
Before running the full MOUSE research code locally, ensure that the following packages are installed:
- [OpenMC](https://github.com/openmc-dev/openmc)
- [WATTS](https://github.com/watts-dev/watts)

## Running Examples Locally

Users can specify reactor design inputs and/or economic parameters for the LTMR, GCMR, and HPMR examples in:
- `examples/watts_exec_LTMR.py`
- `examples/watts_exec_GCMR_Design_A.py`
- `examples/watts_exec_HPMR.py`

A complete detailed bottom-up cost estimation is obtained by running commands such as:
```
python -m examples.watts_exec_LTMR
```
```
python -m examples.watts_exec_GCMR_Design_A
```
```
python -m examples.watts_exec_HPMR
```
Examples of the results are [here](./assets/Ref_Results)

## Citation

If you use MOUSE in technical work, please cite the relevant reports and publications below.

## Relevant Publications

### Foundational Reports and Papers

**[Technoeconomic Evaluation of Microreactor Using Detailed Bottom-up Estimate (Rev. 1)](https://www.osti.gov/biblio/2447366)**
Hanna et al., INL Technical Report, 2024. Develops the transparent bottom-up cost methodology using the MARVEL microreactor as reference — the direct precursor to MOUSE's economic framework.

---

**[Cost Breakdown and Evolution of the MARVEL Microreactor Project](https://www.tandfonline.com/doi/full/10.1080/00295450.2026.2659431)**
Hanna, Strain, Schwartz, and Abou-Jaoude, Nuclear Technology, 2026. Presents a detailed MARVEL microreactor cost breakdown and cost evolution that supports public microreactor cost benchmarking and MOUSE's bottom-up economic framework.

---

**[A Bottom-Up Cost Estimation Tool for Nuclear Microreactors](https://www.osti.gov/biblio/2588465)**
Hanna et al., INL S&T Accomplishment Report, 2025. Describes the MOUSE tool, its architecture, and demonstrates FOAK/NOAK cost estimates and parametric studies for the LTMR, GCMR, and HPMR designs.

---

**[Open-Source Microreactor Design Models for Technoeconomic Assessments](https://www.sciencedirect.com/science/article/pii/S0029549325003875)**
Al-Dawood, Hanna et al., Nuclear Engineering and Design, 2025. Documents the cost correlations, design parameters, and bottom-up assumptions for the LTMR and GCMR models embedded in MOUSE.

---

### Applications and Extensions

**[Techno-Economic Optimization of a Heat-Pipe Microreactor, Part I: Theory and Cost Optimization](https://arxiv.org/pdf/2512.16032)**
Seurin, Price, Nunez — INL/MIT, 2025. Couples MOUSE's LCOE estimation with surrogate modeling and reinforcement learning to optimize HPMR geometry under physics constraints.

---

**[Techno-Economic Optimization of a Heat-Pipe Microreactor, Part II: Multi-Objective Optimization Analysis](https://arxiv.org/pdf/2601.20079)**
Seurin, Price — INL/MIT, 2026. Extends Part I to multi-objective optimization of LCOE and rod-integrated peaking factor using the PEARL algorithm, with MOUSE as the cost engine.

## License

MOUSE is released under the MIT License. See [LICENSE](./LICENSE).
