# Dynamic Systems Simulation Coursework

This repository contains the MATLAB implementations, detailed analytical solutions, simulation results, and technical reports prepared for the graduate course **Dynamic Systems Simulation**, taught by **Dr. Ettefagh**.

The exercises cover the mathematical modeling and numerical simulation of discrete, nonlinear, continuous, and randomly excited mechanical systems. The governing equations are primarily derived using energy methods and Lagrange’s equations and are subsequently implemented and analyzed in MATLAB.

## Repository Contents

Each exercise is organized in a separate folder and may include:

- MATLAB source codes
- Detailed mathematical derivations
- PDF technical reports
- Time-domain simulation results
- Frequency-domain analyses
- Phase-plane diagrams
- FFT and power spectral density results
- A separate README describing the corresponding exercise

## Exercises

### Exercise 1 — Rod–Spring–Damper System with Rotating Unbalance

Modeling and vibration analysis of a two-degree-of-freedom rigid rod supported by springs and a viscous damper and excited by a rotating unbalanced mass.

The exercise includes:

- Derivation of the equations of motion using Lagrange’s equations
- Translational and rotational dynamics
- State-space formulation
- Numerical solution using MATLAB
- Comparison of symmetric and asymmetric spring configurations
- Natural-frequency and resonance analysis
- Time-varying mass and stiffness parameters
- Sensitivity analysis with respect to spring stiffness and damping

### Exercise 2 — Nonlinear Two-DOF Cart–Rod System

Modeling and simulation of a nonlinear two-degree-of-freedom system consisting of a translating mass and a suspended rigid rod with linear and torsional springs and dampers.

The exercise includes:

- Nonlinear Lagrangian modeling
- Translational–rotational dynamic coupling
- Free-vibration analysis
- Impact and pulse excitation
- Harmonic excitation
- Phase-plane analysis
- Investigation of nonlinear softening and response asymmetry
- Identification of translational and torsional resonance regions

### Exercise 3 — Euler–Bernoulli Beam with Concentrated Spring and Damper

Analytical and numerical vibration analysis of a simply supported Euler–Bernoulli beam equipped with a concentrated spring and viscous damper.

The exercise includes:

- Derivation of the governing partial differential equation
- Modal approximation and model reduction
- Construction of modal mass, stiffness, and damping matrices
- Natural-frequency calculation
- Modal convergence analysis
- Free and forced vibration simulations
- Pulse, harmonic, multi-sine, and random excitations
- Fast Fourier Transform analysis
- Investigation of the influence of force and measurement locations

### Exercise 4 — Dynamic Response under Gaussian White-Noise Excitation

Time- and frequency-domain analysis of several mechanical systems subjected to Gaussian white-noise excitation.

The investigated systems include:

1. A rigid rod supported by springs and a damper
2. A nonlinear cart–rod system
3. An Euler–Bernoulli beam with a concentrated spring and damper

The exercise includes:

- Generation of Gaussian white-noise inputs
- Numerical simulation of random vibration responses
- RMS and peak-response evaluation
- Fast Fourier Transform analysis
- Power spectral density estimation using Welch’s method
- Identification of dominant frequencies
- Comparison of spectral peaks with the natural frequencies of each system

## Repository Structure

```text
dynamic-systems-simulation-coursework/
├── Exercise-01-Rod-Spring-Damper-System/
├── Exercise-02-Nonlinear-Cart-Rod-System/
├── Exercise-03-Euler-Bernoulli-Beam/
├── Exercise-04-Random-Vibration-Analysis/
└── README.md
