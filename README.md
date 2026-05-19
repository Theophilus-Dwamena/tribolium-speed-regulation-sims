# Tribolium Speed Regulation Simulations

MATLAB simulations

## Overview
This repository contains three MATLAB simulations that model gap gene expression dynamics along the anterior-posterior (AP) axis of Tribolium castaneum embryos under the speed regulation framework.

## Simulations

### sweeping_pulse_spacetime.m
Simulates a Gaussian transcriptional pulse travelling at constant speed from posterior to anterior. Records intronic and exonic signals at a fixed posterior observation point to show how a spatial wave produces a temporal pulse at a fixed location.

### speed_regulation_three_gradients.m
Simulates a single-gene speed regulation system under three morphogen gradient regimes: gradient-based (static linear), wavefront-based (retracting step), and graded wavefront (smooth retracting step). Produces animated spatial profiles and kymographs for each regime.

### gradient_shape_domain_widths.m
Compares gene expression domain widths produced by a linear versus a concave power-law morphogen gradient under the gradient-based mode. Illustrates how gradient shape affects the spatial uniformity of expression domains along the AP axis.

## Requirements
MATLAB R2020a or later. No additional toolboxes required.

## Usage
Run each .m file independently in MATLAB. 
Each script is self-contained and produces figures directly.
