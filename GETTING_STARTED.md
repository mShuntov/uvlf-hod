# Getting Started with uvlf-hod

This guide will help you get started with the `uvlf-hod` package for modeling high-redshift galaxy populations using the unified model architecture.

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/yourusername/uvlf-hod.git
cd uvlf-hod
```

### Step 2: Install Dependencies

It's recommended to use a virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

Install the package in development mode:

```bash
pip install -e .
```

Or install with development tools:

```bash
pip install -e ".[dev]"
```

### Step 3: Verify Installation

```python
import uvlf_hod
print(uvlf_hod.__version__)
```

## Core Concepts

### The Unified Model Architecture

`uvlf-hod` uses a clean class hierarchy:

```
UVHMRModel (base class)
    ↓ inherits from
HODModel (full model)
```

- **`UVHMRModel`**: Handles UV-halo mass relations (SFR, luminosity)
- **`HODModel`**: Extends UVHMRModel with occupation distributions (LF, bias)

**Key insight**: `HODModel` includes *everything* from `UVHMRModel` plus HOD functionality!

## Quick Tour

### 1. Creating Your First Model

```python
from uvlf_hod import HODModel

# Create a complete model
model = HODModel(
    z=6.0,           # Redshift
    eps0=0.1,        # Star formation efficiency
    Mc=10**11.5,     # Characteristic halo mass [M_sun]
    a=0.6,           # Low-mass slope
    b=0.35,          # High-mass slope
    sigma_UV=0.35,   # UV magnitude scatter [mag]
    Mcut=10**10,     # Satellite cutoff mass [M_sun]
    Msat=10**12.5,   # Satellite normalization [M_sun]
    asat=1.0         # Satellite power-law slope
)

print(model)
```

### 2. Computing UV Luminosity Functions

```python
import numpy as np
import matplotlib.pyplot as plt

# Define magnitude range
MUV = np.linspace(-22, -16, 25)

# Compute luminosity function
phi = model.luminosity_function(MUV)

# Plot
plt.semilogy(MUV, phi)
plt.xlabel('$M_{\\mathrm{UV}}$')
plt.ylabel('$\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.title(f'UV Luminosity Function at z={model.z}')
plt.grid(True, alpha=0.3)
plt.show()
```

### 3. Understanding the Parameters

#### UVHMR Parameters (control star formation)

```python
# eps0: Star formation efficiency peak
# Higher eps0 → more star formation → brighter galaxies
model.update_parameters(eps0=0.15)

# Mc: Characteristic halo mass where efficiency peaks
# Typical: 10^11 - 10^12 M_sun
model.update_parameters(Mc=10**12)

# a, b: Shape of efficiency curve
# a controls low-mass slope, b controls high-mass slope
model.update_parameters(a=0.7, b=0.3)
```

#### HOD Parameters (control occupation)

```python
# sigma_UV: Scatter in M_UV at fixed halo mass
# Larger scatter → more spread in the relation
model.update_parameters(sigma_UV=0.4)

# Mcut: Minimum mass for satellite galaxies
# Only halos with M > Mcut can host satellites
model.update_parameters(Mcut=10**10.5)

# Msat: Normalization for satellite occupation
# Controls how many satellites in massive halos
model.update_parameters(Msat=10**12)

# asat: Power-law slope for satellites
# Steeper slope → more satellites in massive halos
model.update_parameters(asat=1.2)
```

### 4. Using UVHMR Methods

The `HODModel` inherits all methods from `UVHMRModel`:

```python
# Star formation rate for a halo
Mh = 1e11  # M_sun
sfr = model.sfr(Mh)
print(f"SFR = {sfr:.2f} M_sun/yr")

# UV magnitude for a halo
MUV = model.MUV(Mh)
print(f"M_UV = {MUV:.2f}")

# Inverse: halo mass for a UV magnitude
Mh_recovered = model.Mhalo(MUV)
print(f"Recovered: {Mh_recovered:.2e} M_sun")

# Star formation efficiency
epsilon = model.star_formation_efficiency(Mh)
print(f"ε = {epsilon:.3f}")
```

### 5. Using HOD Methods

```python
# Occupation numbers
MUV_thresh = -18
Mh_array = np.logspace(9, 13, 100)

N_cen = model.Ncen(Mh_array, MUV_thresh)  # Central galaxies
N_sat = model.Nsat(Mh_array, MUV_thresh)  # Satellite galaxies
N_tot = model.Ngal(Mh_array, MUV_thresh)  # Total galaxies

# Plot occupation functions
plt.figure(figsize=(8, 6))
plt.plot(np.log10(Mh_array), N_cen, label='Central')
plt.plot(np.log10(Mh_array), N_sat, label='Satellite')
plt.plot(np.log10(Mh_array), N_tot, '--', label='Total')
plt.xlabel('$\\log_{10}(M_h / M_\\odot)$')
plt.ylabel('$\\langle N \\rangle$')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### 6. Galaxy Bias

```python
# Compute bias as a function of magnitude
bias = model.galaxy_bias(MUV)

plt.plot(MUV, bias, linewidth=2)
plt.xlabel('$M_{\\mathrm{UV}}$')
plt.ylabel('Galaxy Bias $b_g$')
plt.title('Galaxy Clustering Bias')
plt.grid(True, alpha=0.3)
plt.show()
```

### 7. Mean Properties

```python
# Properties for galaxies brighter than a threshold
MUV_thresh = -20

mean_mass = model.mean_halo_mass(MUV_thresh)
mean_bias = model.mean_bias(MUV_thresh)

print(f"Galaxies brighter than M_UV = {MUV_thresh}:")
print(f"  Mean halo mass: {10**mean_mass:.2e} M_sun")
print(f"  Mean bias: {mean_bias:.2f}")
```

### 8. Updating Parameters

```python
# Update one or more parameters
model.update_parameters(z=7.0, eps0=0.15)

# The model is ready to use with new parameters
phi_new = model.luminosity_function(MUV)

# Reset if needed
model.update_parameters(z=6.0, eps0=0.1)
```

## Common Workflows

### Workflow 1: Compare Model to Observations

```python
import numpy as np
import matplotlib.pyplot as plt
from uvlf_hod import HODModel

# Observational data (example)
MUV_obs = np.array([-21.5, -20.5, -19.5, -18.5])
phi_obs = np.array([1e-5, 5e-5, 2e-4, 8e-4])
phi_err = phi_obs * 0.3

# Create model
model = HODModel(
    z=6.0, eps0=0.1, Mc=10**11.5, a=0.6, b=0.35,
    sigma_UV=0.35, Mcut=10**10, Msat=10**12.5, asat=1.0
)

# Compute model predictions
MUV_model = np.linspace(-22, -17, 30)
phi_model = model.luminosity_function(MUV_model)

# Plot comparison
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(MUV_obs, phi_obs, yerr=phi_err, fmt='o', 
            label='Observations', capsize=5, markersize=8)
ax.semilogy(MUV_model, phi_model, '-', linewidth=2, label='Model')
ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

### Workflow 2: Parameter Exploration

```python
# Explore effect of varying eps0
eps0_values = [0.05, 0.1, 0.15, 0.2]
colors = plt.cm.viridis(np.linspace(0, 1, len(eps0_values)))

fig, ax = plt.subplots(figsize=(10, 6))

for eps0, color in zip(eps0_values, colors):
    model_var = HODModel(
        z=6.0, eps0=eps0, Mc=10**11.5, a=0.6, b=0.35,
        sigma_UV=0.35, Mcut=10**10, Msat=10**12.5, asat=1.0
    )
    phi = model_var.luminosity_function(MUV_model)
    ax.semilogy(MUV_model, phi, linewidth=2.5, color=color,
                label=f'$\\epsilon_0={eps0}$')

ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

### Workflow 3: Multi-Redshift Analysis

```python
redshifts = [5, 6, 7, 8]
colors = plt.cm.plasma(np.linspace(0, 0.8, len(redshifts)))

fig, ax = plt.subplots(figsize=(10, 6))

# Create one model and update redshift
model = HODModel(
    z=5.0, eps0=0.1, Mc=10**11.5, a=0.6, b=0.35,
    sigma_UV=0.35, Mcut=10**10, Msat=10**12.5, asat=1.0
)

for z_val, color in zip(redshifts, colors):
    model.update_parameters(z=z_val)
    phi = model.luminosity_function(MUV_model)
    ax.semilogy(MUV_model, phi, linewidth=2.5, color=color,
                label=f'z={z_val}')

ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

### Workflow 4: Comprehensive Analysis

```python
# Create model
model = HODModel(
    z=6.0, eps0=0.1, Mc=10**11.5, a=0.6, b=0.35,
    sigma_UV=0.35, Mcut=10**10, Msat=10**12.5, asat=1.0
)

# Compute all quantities
MUV = np.linspace(-22, -16, 25)
phi = model.luminosity_function(MUV)
bias = model.galaxy_bias(MUV)

# Mean properties
mean_mass_bright = model.mean_halo_mass(-20)
mean_mass_faint = model.mean_halo_mass(-17)

# Create comprehensive plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

# Luminosity function
ax1.semilogy(MUV, phi, 'b-', linewidth=2)
ax1.set_xlabel('$M_{\\mathrm{UV}}$')
ax1.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
ax1.set_title('Luminosity Function')
ax1.grid(True, alpha=0.3)

# Galaxy bias
ax2.plot(MUV, bias, 'r-', linewidth=2)
ax2.set_xlabel('$M_{\\mathrm{UV}}$')
ax2.set_ylabel('$b_g$')
ax2.set_title('Galaxy Bias')
ax2.grid(True, alpha=0.3)

# UVHMR
Mh = np.logspace(9, 13, 100)
MUV_hmr = model.MUV(Mh)
ax3.plot(np.log10(Mh), MUV_hmr, 'g-', linewidth=2)
ax3.set_xlabel('$\\log_{10}(M_h / M_\\odot)$')
ax3.set_ylabel('$M_{\\mathrm{UV}}$')
ax3.invert_yaxis()
ax3.set_title('UV-Halo Mass Relation')
ax3.grid(True, alpha=0.3)

# Occupation
N_cen = model.Ncen(Mh, -18)
N_sat = model.Nsat(Mh, -18)
ax4.plot(np.log10(Mh), N_cen, label='Central', linewidth=2)
ax4.plot(np.log10(Mh), N_sat, label='Satellite', linewidth=2)
ax4.set_xlabel('$\\log_{10}(M_h / M_\\odot)$')
ax4.set_ylabel('$\\langle N \\rangle$')
ax4.set_title('Occupation (M_UV < -18)')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

print("\nResults Summary:")
print(f"  Bright galaxies (M_UV < -20): <M_h> = {10**mean_mass_bright:.2e} M_sun")
print(f"  Faint galaxies (M_UV < -17): <M_h> = {10**mean_mass_faint:.2e} M_sun")
```

## Advanced Topics

### Using Only UVHMR (No HOD)

If you only need UV-halo mass relations without occupation distributions:

```python
from uvlf_hod import UVHMRModel

uvhmr = UVHMRModel(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)

# UVHMR methods available
MUV = uvhmr.MUV(1e11)
sfr = uvhmr.sfr(1e11)

# HOD methods NOT available (raises error)
# uvhmr.luminosity_function(MUV)  # ❌ AttributeError
```

### Custom Cosmology

```python
from uvlf_hod.config import CosmologyConfig

# Create custom cosmology
custom_cosmo = CosmologyConfig(H0=67.0, Om0=0.32, Ob0=0.05)

# Models will use this cosmology
```

### Redshift Evolution of Parameters

```python
from uvlf_hod.models.parametrization import eps0_fz

# Define evolution
z_array = np.linspace(4, 10, 20)
eps0_values = eps0_fz(z_array, deps_dz=0.01, eps_off=0.05)

# Create models at each redshift
models = [
    HODModel(z=z, eps0=eps0, Mc=10**11.5, a=0.6, b=0.35,
             sigma_UV=0.35, Mcut=10**10, Msat=10**12.5, asat=1.0)
    for z, eps0 in zip(z_array, eps0_values)
]
```

## Troubleshooting

### Issue: Import errors

**Solution**: Make sure all dependencies are installed:
```bash
pip install -r requirements.txt
```

### Issue: Slow calculations

**Solution**: The first calculation creates the HMF cache. Subsequent calculations are much faster. Reduce `num_points` in `get_halo_mass_function` for faster (but less accurate) results.

### Issue: NaN values in results

**Solution**: Check parameter ranges. Very extreme values can cause numerical issues. Ensure masses are reasonable (10^9 - 10^14 M_sun).

## Next Steps

- Read the full [API documentation](docs/api/)
- Check out the [Jupyter notebook](examples/basic_usage.ipynb)
- Explore [additional examples](examples/)
- Read the [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) if upgrading from old version

## Getting Help

- Open an issue on [GitHub](https://github.com/yourusername/uvlf-hod/issues)
- Check the examples directory
- Contact: your.email@example.com