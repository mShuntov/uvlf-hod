# Getting Started with halogal

This guide will help you get started with the `halogal` package for modeling high-redshift galaxy populations using the unified model architecture.

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/mshuntov/halogal.git
cd halogal
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
import halogal
print(halogal.__version__)
```

## Core Concepts

### The Unified Model Architecture

`halogal` uses a clean class hierarchy:

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
from halogal import HODModel

# Create a complete model using fitted parameters from Shuntov+2025
model = HODModel(
    z=6.0,           # Redshift
    eps0=0.19,       # Star formation efficiency
    Mc=10**11.64,    # Characteristic halo mass [M_sun]
    a=0.69,          # Low-mass slope (beta)
    b=0.65,          # High-mass slope (gamma)
    sigma_UV=0.69,   # UV magnitude scatter [mag]
    Mcut=10**9.57,   # Satellite cutoff mass [M_sun]
    Msat=10**12.65,  # Satellite normalization [M_sun]
    asat=0.85        # Satellite power-law slope
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

The default parameters are taken at z~5.4 from Shuntov+2025.

#### UVHMR Parameters (control star formation)

```python
# eps0: Star formation efficiency peak (default: 0.19)
# Higher eps0 → more star formation → brighter galaxies

# Mc: Characteristic halo mass where efficiency peaks (default: 10^11.64 M_sun)

# a, b: Shape of efficiency curve (defaults: a=0.69, b=0.65)
# a controls low-mass slope, b controls high-mass slope
```

#### HOD Parameters (control occupation)

```python
# sigma_UV: Scatter in M_UV at fixed halo mass (default: 0.69 mag)

# Mcut: Minimum mass for satellite galaxies (default: 10^9.57 M_sun)
# Only halos with M > Mcut can host satellites

# Msat: Normalization for satellite occupation (default: 10^12.65 M_sun)
# Controls how many satellites in massive halos

# asat: Power-law slope for satellites (default: 0.85)
```

#### Using Defaults from Config

```python
from halogal.config import DEFAULT_HOD_PARAMS

# All defaults from Shuntov+2025
print(f"eps0 = {DEFAULT_HOD_PARAMS['eps0']}")           # 0.19
print(f"log(Mc) = {DEFAULT_HOD_PARAMS['logMc']}")       # 11.64
print(f"a = {DEFAULT_HOD_PARAMS['a']}")                 # 0.69
print(f"b = {DEFAULT_HOD_PARAMS['b']}")                 # 0.65
print(f"sigma_UV = {DEFAULT_HOD_PARAMS['sigma_UV']}")   # 0.69
print(f"log(Mcut) = {DEFAULT_HOD_PARAMS['Mcut']}")      # 9.57
print(f"log(Msat) = {DEFAULT_HOD_PARAMS['Msat']}")      # 12.65
print(f"asat = {DEFAULT_HOD_PARAMS['asat']}")           # 0.85
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
```

### 5. Using HOD Methods

```python
# Occupation numbers
MUV_thresh = -18
Mh_array = np.logspace(9, 13, 100)

N_cen = model.Ncen(Mh_array, MUV_thresh)  # Central galaxies
N_sat = model.Nsat(Mh_array, MUV_thresh)  # Satellite galaxies
N_tot = model.Ngal(Mh_array, MUV_thresh)  # Total galaxies

# Plot
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

## Common Workflows

### Workflow 1: Compare Model to Observations

```python
from bouwens21_data import bouwens21, redshift_centers

# Load Bouwens+2021 data
obs = bouwens21['z6']
z_obs = redshift_centers['z6']

# Create model
model = HODModel(
    z=z_obs, eps0=0.19, Mc=10**11.64, a=0.69, b=0.65,
    sigma_UV=0.69, Mcut=10**9.57, Msat=10**12.65, asat=0.85
)

# Compute model
MUV_model = np.linspace(-23, -15, 50)
phi_model = model.luminosity_function(MUV_model)

# Plot comparison
fig, ax = plt.subplots(figsize=(10, 7))
ax.errorbar(obs['M_AB'], obs['Fi_k'], yerr=obs['Fi_k_error'],
            fmt='o', markersize=8, capsize=4, label='Bouwens+2021')
ax.semilogy(MUV_model, phi_model, 'r-', linewidth=2.5, label='Model')
ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.title(f'z = {z_obs}', fontsize=16)
plt.show()
```

### Workflow 2: Redshift Evolution of Parameters

```python
from halogal.models.parametrization import eps0_fz, Mc_fz, a_fz, b_fz
from halogal.config import DEFAULT_REDSHIFT_EVOLUTION

# Define redshift range
z_array = np.linspace(4, 8, 20)

# Get evolved parameters using fitted evolution
eps0_z = eps0_fz(
    z_array,
    deps_dz=DEFAULT_REDSHIFT_EVOLUTION['d_eps0_dz'],  # 0.02
    eps_off=DEFAULT_REDSHIFT_EVOLUTION['C_eps0']      # 0.09
)

Mc_z = 10**Mc_fz(
    z_array,
    dMc_dz=DEFAULT_REDSHIFT_EVOLUTION['d_logMc_dz'],  # 0.14
    Mc_off=DEFAULT_REDSHIFT_EVOLUTION['C_logMc']       # 10.86
)

a_z = a_fz(
    z_array,
    da_dz=DEFAULT_REDSHIFT_EVOLUTION['d_a_dz'],        # 0.01
    a_off=DEFAULT_REDSHIFT_EVOLUTION['C_a']            # 0.62
)

b_z = b_fz(
    z_array,
    db_dz=DEFAULT_REDSHIFT_EVOLUTION['d_b_dz'],        # -0.06
    b_off=DEFAULT_REDSHIFT_EVOLUTION['C_b']            # 0.97
)

# Plot parameter evolution
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

axes[0, 0].plot(z_array, eps0_z, 'b-', linewidth=2)
axes[0, 0].set_ylabel('$\\epsilon_0$', fontsize=14)
axes[0, 0].set_xlabel('Redshift', fontsize=12)
axes[0, 0].grid(True, alpha=0.3)

axes[0, 1].plot(z_array, np.log10(Mc_z), 'r-', linewidth=2)
axes[0, 1].set_ylabel('$\\log_{10}(M_c / M_\\odot)$', fontsize=14)
axes[0, 1].set_xlabel('Redshift', fontsize=12)
axes[0, 1].grid(True, alpha=0.3)

axes[1, 0].plot(z_array, a_z, 'g-', linewidth=2)
axes[1, 0].set_ylabel('$a$ (low-mass slope)', fontsize=14)
axes[1, 0].set_xlabel('Redshift', fontsize=12)
axes[1, 0].grid(True, alpha=0.3)

axes[1, 1].plot(z_array, b_z, 'm-', linewidth=2)
axes[1, 1].set_ylabel('$b$ (high-mass slope)', fontsize=14)
axes[1, 1].set_xlabel('Redshift', fontsize=12)
axes[1, 1].grid(True, alpha=0.3)

plt.suptitle('Parameter Evolution with Redshift', fontsize=16)
plt.tight_layout()
plt.show()

# Create models with evolved parameters
colors = plt.cm.viridis(np.linspace(0, 1, len(z_array)))
fig, ax = plt.subplots(figsize=(10, 7))

for i, (z, eps0, Mc, a, b, color) in enumerate(zip(z_array, eps0_z, Mc_z, a_z, b_z, colors)):
    if i % 4 == 0:  # Plot every 4th redshift
        model = HODModel(
            z=z, eps0=eps0, Mc=Mc, a=a, b=b,
            sigma_UV=0.69, Mcut=10**9.57, Msat=10**12.65, asat=0.85
        )
        MUV = np.linspace(-22, -16, 30)
        phi = model.luminosity_function(MUV)
        ax.semilogy(MUV, phi, color=color, linewidth=2, label=f'z={z:.1f}')

ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
plt.title('Luminosity Function Evolution', fontsize=16)
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
    z=5.0, eps0=0.19, Mc=10**11.64, a=0.69, b=0.65,
    sigma_UV=0.69, Mcut=10**9.57, Msat=10**12.65, asat=0.85
)

MUV_model = np.linspace(-22, -16, 30)

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


### Custom Cosmology

```python
from halogal.config import CosmologyConfig

# Create custom cosmology
custom_cosmo = CosmologyConfig(H0=67.0, Om0=0.32, Ob0=0.05)
```

## Troubleshooting

### Issue: Import errors

**Solution**: Make sure all dependencies are installed:
```bash
pip install -r requirements.txt
```

### Issue: Slow calculations

**Solution**: The first calculation creates the HMF cache. Subsequent calculations are much faster.

### Issue: NaN values in results

**Solution**: Check parameter ranges. Very extreme values can cause numerical issues.

## Next Steps

- Read the full [API documentation](docs/api/)
- Check out the [Jupyter notebook](examples/basic_usage.ipynb)
- Explore [additional examples](examples/)

## Getting Help

- Open an issue on [GitHub](https://github.com/mshuntov/halogal/issues)
- Check the examples directory
- Read the paper: Shuntov et al. 2025, A&A 699 A231
- Contact: marko.shuntov@nbi.ku.dk