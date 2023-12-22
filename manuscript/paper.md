---
title: 'amisrsynthdata: A Python package for generating synthetic data for the Advanced Modular Incoherent Scatter Radars'
tags:
  - Python
  - AMISR
  - ionosphere
authors:
  - name: Leslie J. Lamarche
    orcid: 0000-0001-7098-0524
    affiliation: 1
affiliations:
  - name: SRI International, Menlo Park, CA, USA
    index: 1
date: 15 December 2023
bibliography: paper.bib
---

# Summary

Advanced Modular Incoherent Scatter Radars (AMISR) are a class of ground-based radar facility that are used to study plasma dynamics in the ionosphere [@Kelly2009].  These remotely measure the electron density, $N_e$, ion temperature, $T_i$, electron temperature, $T_e$, and line-of-site component of the drift velocity, $V_{los}$, of the plasma in the ionosphere.  AMISR utilizes phased-array technology to steer the radar beam electronically, which allows it to change look directions many times per second.  This allows the radar to both effectively collect data from different beam positions simultaneously and have a great degree of flexibility in terms of beam positions and mode designs.  Most standard modes use between 4 and 50 beams, roughly evenly spread across the radar field-of-view.  Some specialized modes use tight clusters of beams, or other unusual configurations to observe specific phenomena. This flexibility, as well as the radars' locations in regions that are key for magnetosphere-ionosphere coupling, make them a heavily utilized and invaluable resource for space physics research.

The `amisrsynthdata` package produces synthetic data for AMISR.  It is written purely in python, with use of common numeric and array manipulation libraries [@numpy; @h5py; @pymap3d; @matplotlib; @cartopy].  Users create a configuration file containing information both about the radar operational mode and the ionospheric states, which the package then uses to create a synthetic data file showing what the data would be expected to look like if the radar were to measure the specified phenomena in the specified mode.  This is useful for several things related to ensuring efficient and effective use of the facilities.

1. **Validation of high-level data products:** AMISR data are often used in sophisticated inversion and interpolation procedures to generate high level data products.  These include things such as generating vector velocities from the line-of-site measurements [@Heinselman2008; @Nicolls2014], 3D volumetric interpolations [@Lamarche2020], and inversion of the precipitating energetic particle spectra [@Semeter2005].  When developing these algorithms, it is important to have "truth" data to validate against, which is only possible with synthetic data.

2. **Designing and optimizing observational modes:**  The `amisrsynthdata` package lets users experiment with different beam positions and modes, which allows you to determine the optimal configuration for observing particular features or phenomena.  This is particularly useful when coordinating with rare events such as rocket or satellite conjunctions and there will not be the opportunity to iterate over the mode through normal observations.

3. **Determining if a particular phenomena is observable:** Sometime it is unclear if an ionospheric structure can be seen in the radar data, or what it is expected to look like.  This is relevant to determine if a radar mode under-resolves a structure, which can lead to aliasing and bias interpretation if not accounted for.  This is also useful if you are looking for a particular kind of event in the collected database, but need to understand the expected signatures of it in the radar data.

Although there are only a small number of people who routinely design and optimize radar operation modes, this work is done in service to the broader ionospheric science research community.  Furthermore, the user community routinely experiments with new ways to utilize the radar data and extract more information from it.  AMISR is a widely utilized facility that continues to contribute significantly to space physics research [@decadalsurvey2013].

![Example of summary figure output of `amisrsynthdata`.  This shows electron density for a traveling ionospheric disturbance or large-scale wave propagating in the F-region.  The top row of panels show the "truth" density at different altitude slices as well as the beam locations in each.  The bottom panel shows a "range-time-intensity" plot of synthetic measured electron density over time in a specific beam.  The right panel shows a 3D view of synthetic measurements in all beams at a particular time.  The beam used in the bottom panel is indicated in pink in the top panels and the time used in the top and right panels is indicated by the pink line in the bottom panel.](../docs/synthdata_summary_ne.png)

# Acknowledgements

AMISR PI Dr. Asti Bhatt and former PI Dr. Roger Varney are graciously thanked for their support of this work.  The AMISR facilities are funded by the National Science Foundation through cooperative agreement AGS-1840962 to SRI International. The development of `amisrsynthdata` was funded in part by NASA awards 80NSSC21K0458, 80NSSC21K1354, and 80NSSC21K1318 and NSF award 2027300.

# References
