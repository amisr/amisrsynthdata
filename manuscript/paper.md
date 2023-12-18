---
title: 'amisrsynthdata: A Python package for generating synthetic data for the Advanced Modular Incoherent Scatter Radars'
tags:
  - Python
  - AMISR
  - ionosphere
authors:
  - name: Leslie J. Lamarche
    orcid: 0000-1234
    affiliation: 1
  - name: Asti Bhatt
    orcid: 0000-1234
    affiliation: 1
affiliations:
  - name: SRI International, Menlo Park, CA, USA
  - index: 1
date: 15 December 2023
bibliography: paper.bib
---

# Summary

Advanced Modular Incoherent Scatter Radars (AMISR) are a class of ground-based radar facility that are used to study plasma dynamics in the ionosphere.  These remotely measure the electron density, $N_e$, ion temperature, $T_i$, electron temperature, $T_e$, and line-of-site component of the drift velocity, $V_{los}$, of the plasma in the ionosphere.  AMISR utilizes phased-array technology to be able to stear the radar beam electronically, which allows it to change look directions muliple times per second.  This allows the radar to both collect data from different beam positions effectively simultaniously and have a great degree of flexibility in terms of beam positions and mode designs.  For instance, some common modes include WorldDay, with 11 beams roughly evently spread, or Themis, with 23 beams roughly evenly spread, but thre are also more unusual modes used for specialized observations, such as the GPS mode which includes 5 beams in a very tight cross to make detailed measurements of gradients along a GNSS line-of-sight through the ionosphere. This flexibility, as well as the radar's positioning in regions that are key for magnetosphere-ionosphere coupling, make them a heavily utilized and invaluable resource for splace physics research.

The `amisrsynthdata` package produces synthetic data for AMISR.  It is written purely in python, with use of common numeric and array manipulation libraries.  Users create a configuration file, which contains information both about the radar operational mode and the ionospheric states, which the package then uses to create a synthetic data file showing what the data would be expected to look like if the radar were to measure the specified phenomena in the specified mode.  This is useful for several things related to ensuring efficient and effective opperations of the radar.

1. Validation of high-level data products. AMISR data are often used in sophisticated inversion and interpolation procedures to generate high level data products.  These include things such as generating 3D vector velocities from the line-of-site vectors measured by the radar and 3D volumetric interpolations.  When developing these algorithms, it is important to have "truth" data to validate against, which is only possible with synthetic data.

2. Designing and optimizing observational modes for specific purposes.  The `amisrsynthdata` package lets users experiment with different beam positions and configurations, which allows you to determine the optimal configuration for a particular feature.  This is particularly useful when the event in question is rare and there will not be the opportunity to iterate over the mode through normal observations.

3. Determine if a particular feaure is observable, and if so, what it is expected to look like in the radar data.  This is useful if you are looking for a particular kind of event in the collected database, but need to understand the signatures of it in the radar data.

Although there are only a small number of people who routinely design and optimize radar operation modes, this work is done in service to the broader ionospheric science research community.  Furthermore, the user community routinely experiments with new ways to utilize the radar data and extract more information from it, which are alternative uses of the `amisrsynthdata` package.

# Acknowledgements

The develpment of `amisrsynthdata` was funded as part of several research grants, including...  The AMISR failities are operated by SRI International under cooperative agreement ... from NSF.

# References
