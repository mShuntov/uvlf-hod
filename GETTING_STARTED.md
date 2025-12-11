# Getting Started with uvlf-hod

This guide will help you get started with the `uvlf-hod` package for modeling high-redshift galaxy populations.

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/mshuntov/uvlf-hod.git
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

## Quick Tour

### 1. Computing a UV Luminosity Function

The main use case is computing UV luminosity functions at high redshift:

```python
import numpy as np
import matplotlib.pyplot as plt
from uvlf_hod import compute_luminosity_function

# Set up parameters
z = 6.0
MUV = np.linspace(-22, -16, 25)

# UVHMR parameters: [eps0, log10(Mc), a, b]
uvhmr_params = [0.1, 11.5, 0.6, 0.35]

# Satellite parameters: [log10(Mcut), log10(Msat), asat]
nsat_params = [10.0, 12.5, 1.0]

# UV scatter
sigma_UV = 0.35

# Compute
phi = compute_luminosity_function(MUV, z, uvhmr_params, nsat_params, sigma_UV)

# Plot
plt.semilogy(MUV, phi)
plt.xlabel('$M_{UV}$')
plt.ylabel('$\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
plt.show()
```

### 2. Understanding the Parameters

#### UVHMR Parameters

- **eps0**: Star formation efficiency normalization (typical: 0.05-0.2)
- **log10(Mc)**: Log of characteristic halo mass (typical: 11-12)
- **a**: Low-mass slope (typical: 0.4-0.8)
- **b**: High-mass slope (typical: 0.2-0.5)

#### HOD Parameters

- **log10(Mcut)**: Log cutoff mass for satellites (typical: 9-11)
- **log10(Msat)**: Log satellite normalization mass (typical: 11-13)
- **asat**: Satellite power-law slope (typical: 0.8-1.2)

#### Other Parameters

- **sigma_UV**: UV magnitude scatter (typical: 0.2-0.5 mag)
- **add_dust**: Whether to include dust attenuation (default: True)

### 3. UV-Halo Mass Relation

Explore the relationship between halo mass and UV magnitude:

```python
from uvlf_hod.core import UVHaloMassRelation

# Create relation
uvhmr = UVHaloMassRelation(z=6.0, eps0=0.1, Mc=1e12, a=0.6, b=0.35)

# Forward: Halo mass to UV magnitude
Mh = 1e11  # Solar masses
MUV = uvhmr.MUV(Mh)
print(f"M_h = {Mh:.2e} -> M_UV = {MUV:.2f}")

# Backward: UV magnitude to halo mass
Mh_recovered = uvhmr.Mhalo(MUV)
print(f"M_UV = {MUV:.2f} -> M_h = {Mh_recovered:.2e}")

# Get SFR
sfr = uvhmr.sfr(Mh)
print(f"SFR = {sfr:.1f} M_sun/yr")
```

### 4. Galaxy Bias

Calculate galaxy clustering bias:

```python
from uvlf_hod import compute_galaxy_bias

bias = compute_galaxy_bias(MUV, z, uvhmr_params, nsat_params, sigma_UV)

plt.plot(MUV, bias)
plt.xlabel('$M_{UV}$')
plt.ylabel('Galaxy Bias')
plt.show()
```

### 5. Working with Cosmology

Customize the cosmology:

```python
from uvlf_hod.config import CosmologyConfig

# Create custom cosmology
cosmo = CosmologyConfig(H0=67.0, Om0=0.32, Ob0=0.05)

# Use with halo mass function
from uvlf_hod.cosmology import HaloMassFunction

hmf = HaloMassFunction(z=6.0, M_min=8, M_max=14)
print(f"Mass range: {hmf.log10_mass.min():.1f} to {hmf.log10_mass.max():.1f}")
```

### 6. Redshift Evolution

Model parameter evolution with redshift:

```python
from uvlf_hod.models.parametrization import eps0_fz, Mc_fz

# Define evolution
z_array = np.linspace(4, 10, 20)

# Evolve eps0: eps0(z) = deps_dz * z + eps_off
eps0_evolution = eps0_fz(z_array, deps_dz=0.01, eps_off=0.05)

# Evolve Mc: Mc(z) = dMc_dz * z + Mc_off
Mc_evolution = Mc_fz(z_array, dMc_dz=0.1, Mc_off=11.0)

plt.plot(z_array, eps0_evolution, label='$\epsilon_0(z)$')
plt.xlabel('Redshift')
plt.ylabel('$\epsilon_0$')
plt.legend()
plt.show()
```

## Common Workflows

### Workflow 1: Compare Model to Observations

```python
import numpy as np
import matplotlib.pyplot as plt
from uvlf_hod import compute_luminosity_function

# Observational data (example)
MUV_obs = np.array([-21.5, -20.5, -19.5, -18.5])
phi_obs = np.array([1e-5, 5e-5, 2e-4, 8e-4])
phi_err = phi_obs * 0.3

# Model
z = 6.0
MUV_model = np.linspace(-22, -17, 30)
uvhmr_params = [0.1, 11.5, 0.6, 0.35]
nsat_params = [10.0, 12.5, 1.0]
sigma_UV = 0.35

phi_model = compute_luminosity_function(
    MUV_model, z, uvhmr_params, nsat_params, sigma_UV
)

# Plot comparison
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(MUV_obs, phi_obs, yerr=phi_err, fmt='o', 
            label='Observations', capsize=5)
ax.semilogy(MUV_model, phi_model, '-', label='Model')
ax.set_xlabel('$M_{\\mathrm{UV}}$')
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
```

### Workflow 2: Parameter Exploration

```python
# Explore effect of varying eps0
eps0_values = [0.05, 0.1, 0.15, 0.2]

fig, ax = plt.subplots(figsize=(8, 6))

for eps0 in eps0_values:
    uvhmr_params = [eps0, 11.5, 0.6, 0.35]
    phi = compute_luminosity_function(
        MUV_model, z, uvhmr_params, nsat_params, sigma_UV
    )
    ax.semilogy(MUV_model, phi, label=f'$\\epsilon_0={eps0}$')

ax.set_xlabel('$M_{\\mathrm{UV}}$')
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
```

### Workflow 3: Multi-redshift Analysis

```python
redshifts = [5, 6, 7, 8]
colors = plt.cm.viridis(np.linspace(0, 1, len(redshifts)))

fig, ax = plt.subplots(figsize=(8, 6))

for z, color in zip(redshifts, colors):
    phi = compute_luminosity_function(
        MUV_model, z, uvhmr_params, nsat_params, sigma_UV
    )
    ax.semilogy(MUV_model, phi, color=color, label=f'z={z}')

ax.set_xlabel('$M_{\\mathrm{UV}}$')
ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
```

## Advanced Topics

### Custom Halo Mass Functions

```python
from uvlf_hod.cosmology import get_halo_mass_function

# Try different models
models = ['watson13', 'tinker08', 'reed07']

fig, ax = plt.subplots()

for model in models:
    log10_M, hmf = get_halo_mass_function(
        z=6.0, M_min=8, M_max=14, model=model
    )
    ax.semilogy(log10_M, hmf, label=model)

ax.set_xlabel('$\\log_{10}(M/M_\\odot)$')
ax.set_ylabel('$dn/d\\ln M$ [Mpc$^{-3}$]')
ax.legend()
plt.show()
```

### Dust Treatment

```python
from uvlf_hod.luminosity import dust_attenuation, beta_color

# Explore dust attenuation
z = 6.0
MUV_range = np.linspace(-22, -16, 50)

A_UV = dust_attenuation(z, MUV_range)
beta = beta_color(z, MUV_range)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(MUV_range, A_UV)
ax1.set_xlabel('$M_{\\mathrm{UV}}$')
ax1.set_ylabel('$A_{\\mathrm{UV}}$ [mag]')
ax1.set_title('Dust Attenuation')

ax2.plot(MUV_range, beta)
ax2.set_xlabel('$M_{\\mathrm{UV}}$')
ax2.set_ylabel('$\\beta$')
ax2.set_title('UV Slope')

plt.tight_layout()
plt.show()
```

## Troubleshooting

### Issue: Import errors

**Solution**: Make sure all dependencies are installed:
```bash
pip install -r requirements.txt
```

### Issue: Slow calculations

**Solution**: Reduce number of mass bins or use parallel processing:
```python
from uvlf_hod.config import HMF_NUM_POINTS
# Reduce from default 2048 to 1024 for faster (but less accurate) results
```

### Issue: NaN values in results

**Solution**: Check parameter ranges. Very extreme values can cause numerical issues.

## Next Steps

- Read the full [API documentation](docs/api/)
- Check out more [examples](examples/)
- Join the community discussions on GitHub
- Contribute improvements via pull requests

## Getting Help

- Open an issue on [GitHub](https://github.com/yourusername/uvlf-hod/issues)
- Check the [FAQ](docs/faq.md)
- Contact: your.email@example.com
