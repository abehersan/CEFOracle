from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np
from CrystalField import CrystalField
from scipy.integrate import trapezoid


def calc_SM(HC, T):
    SM = np.zeros(len(HC))
    for i in range(len(HC)):
        SM[i] = trapezoid(HC[0:i] / T[0:i], T[0:i], dx=0.5)
    return SM


def CEF_comparison():
    cef = CrystalField('Yb', 'C6v')
    cef["B20"] = 0.5622
    cef["B40"] = 1.6087e-5
    cef["B60"] = 6.412e-7
    cef["B66"] = -8.324e-6
    
    fig, axes = plt.subplots(2, 2, figsize=(11.7, 8.3))
    ax_ins, ax_mag, ax_chi, ax_hc = axes.ravel()
    
    ### INS X-section ###
    cef.PeakShape = 'Gaussian'
    cef.FWHM = 0.5
    EN = np.linspace(-15, 15, 701)
    for T in [10.0, 50.0, 200.0]:
        cef.Temperature = T
        ins = cef.getSpectrum(EN)
        ax_ins.plot(*ins, label=f"T = {T}K", lw=1)
    ax_ins.set_xlim(-15, 15)
    ax_ins.legend()
    ax_ins.set_xlabel("E [meV]")
    ax_ins.set_ylabel("I(Q, E) [arb. units]")
    
    ### Magnetisation ###
    Bs = np.linspace(0, 12, 101)
    Mu_para = cef.getMagneticMoment(Temperature=1.5, Hmag=Bs, Hdir=[0,0,1], Unit="bohr")
    Mu_perp = cef.getMagneticMoment(Temperature=1.5, Hmag=Bs, Hdir=[1,0,0], Unit="bohr")
    ax_mag.plot(*Mu_para, label="B parallel c", lw=1)
    ax_mag.plot(*Mu_perp, label="B perp c", lw=1)
    ax_mag.set_xlim(0, 12)
    ax_mag.set_ylim(0, 3.5)
    ax_mag.legend()
    ax_mag.set_xlabel("B [T]")
    ax_mag.set_ylabel("mu [muB / ion]")
    
    ### Static susceptibility ###
    Ts = np.linspace(0, 300, 701)
    chi_para = cef.getSusceptibility(Ts, Hdir=[0,0,1], Unit="CGS")
    chi_perp = cef.getSusceptibility(Ts, Hdir=[1,0,0], Unit="CGS")
    chi_powd = cef.getSusceptibility(Ts, Hdir="powder", Unit="CGS")
    ax_chi.plot(Ts, 1 / chi_para[1], label="B parallel c", lw=1)
    ax_chi.plot(Ts, 1 / chi_perp[1], label="B perp c", lw=1)
    ax_chi.plot(Ts, 1 / chi_powd[1], label="powder", lw=1)
    ax_chi.set_xlim(0, 300)
    ax_chi.set_ylim(0, 155)
    ax_chi.legend()
    ax_chi.set_xlabel("T [K]")
    ax_chi.set_ylabel("1/chi [emu/cm^3]")
    
    ### Schottky heat capacity and magnetic entropy ###
    Ts = np.linspace(0.5, 300, 701)
    HC = cef.getHeatCapacity(Ts)
    SM = calc_SM(HC[1], Ts)
    ax_hc.plot(Ts, HC[1], label="HC", lw=1)
    ax_hc.plot(Ts, SM, label="SM", lw=1)
    ax_hc.set_xlim(0, 300)
    ax_hc.set_ylim(0, 15)
    ax_hc.legend()
    ax_hc.set_xlabel("T [K]")
    ax_hc.set_ylabel("HC [J/mol/K]")
    
    plt.tight_layout()
    plt.show()
    
    return None
    
    
CEF_comparison()