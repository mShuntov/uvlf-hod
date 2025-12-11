"""
Basic usage example for uvlf-hod package.

This script demonstrates how to:
1. Set up model parameters
2. Compute UV luminosity functions
3. Calculate galaxy bias
4. Plot results
"""

import numpy as np
import matplotlib.pyplot as plt
from uvlf_hod import compute_luminosity_function, compute_galaxy_bias
from uvlf_hod.core import UVHaloMassRelation


def main():
    """Run basic examples."""
    
    # =================================================================
    # Example 1: UV Luminosity Function at z=6
    # =================================================================
    print("Computing UV Luminosity Function at z=6...")
    
    # Define redshift
    z = 6.0
    
    # Define UV magnitude range
    MUV = np.linspace(-22, -16, 25)
    
    # Set up UVHMR parameters [eps0, log10(Mc), a, b]
    uvhmr_params = [0.1, 11.5, 0.6, 0.35]
    
    # Satellite parameters [log10(Mcut), log10(Msat), asat]
    nsat_params = [10.0, 12.5, 1.0]
    
    # UV magnitude scatter
    sigma_UV = 0.35
    
    # Compute luminosity function
    phi = compute_luminosity_function(
        MUV, z, uvhmr_params, nsat_params, sigma_UV
    )
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.semilogy(MUV, phi, 'b-', linewidth=2, label=f'z={z}')
    ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
    ax.set_ylabel('$\\Phi$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=14)
    ax.set_title('UV Luminosity Function', fontsize=16)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('luminosity_function.png', dpi=150)
    print("Saved luminosity_function.png")
    
    # =================================================================
    # Example 2: Galaxy Bias
    # =================================================================
    print("\nComputing galaxy bias...")
    
    # Compute bias for same magnitude bins
    bias = compute_galaxy_bias(MUV, z, uvhmr_params, nsat_params, sigma_UV)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(MUV, bias, 'r-', linewidth=2, label=f'z={z}')
    ax.set_xlabel('$M_{\\mathrm{UV}}$', fontsize=14)
    ax.set_ylabel('Galaxy Bias $b_g$', fontsize=14)
    ax.set_title('Galaxy Clustering Bias', fontsize=16)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('galaxy_bias.png', dpi=150)
    print("Saved galaxy_bias.png")
    
    # =================================================================
    # Example 3: UV-Halo Mass Relation
    # =================================================================
    print("\nExploring UV-Halo Mass Relation...")
    
    # Create UVHMR instance
    eps0, log_Mc, a, b = uvhmr_params
    uvhmr = UVHaloMassRelation(
        z=z, eps0=eps0, Mc=10**log_Mc, a=a, b=b, add_dust=True
    )
    
    # Halo mass range
    Mh_array = np.logspace(9, 13, 100)
    
    # Get UV magnitudes
    MUV_from_Mh = uvhmr.MUV(Mh_array)
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(np.log10(Mh_array), MUV_from_Mh, 'g-', linewidth=2)
    ax.set_xlabel('$\\log_{10}(M_h / M_\\odot)$', fontsize=14)
    ax.set_ylabel('$M_{\\mathrm{UV}}$', fontsize=14)
    ax.set_title('UV-Halo Mass Relation', fontsize=16)
    ax.invert_yaxis()  # Brighter magnitudes on top
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('uvhmr.png', dpi=150)
    print("Saved uvhmr.png")
    
    # =================================================================
    # Example 4: Print some key values
    # =================================================================
    print("\n" + "="*60)
    print("Key Results:")
    print("="*60)
    
    # Find characteristic values
    idx_bright = np.argmin(np.abs(MUV + 20))
    print(f"\nAt M_UV = {MUV[idx_bright]:.1f}:")
    print(f"  Φ = {phi[idx_bright]:.2e} Mpc^-3 mag^-1")
    print(f"  Bias = {bias[idx_bright]:.2f}")
    
    # Halo mass for bright galaxy
    Mh_bright = uvhmr.Mhalo(MUV[idx_bright])
    print(f"  Typical halo mass = {Mh_bright:.2e} M_sun")
    
    # SFR
    sfr_bright = uvhmr.sfr(Mh_bright)
    print(f"  Star formation rate = {sfr_bright} M_sun/yr")
    
    print("\n" + "="*60)
    print("Example completed successfully!")
    print("="*60)


if __name__ == "__main__":
    main()
