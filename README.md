# Polymer models of mitotic chromosomes, v2.0

This repository contains the scripts that generate polymer models of mitotic chromosomes published in the manuscript "Rules of engagement for condensins and cohesins guide mitotic chromosome formation" by Samejima, Gibcus et al (Mirny/Dekker/Goloborodko/Earnshaw labs) ([bioRxiv](https://www.biorxiv.org/content/10.1101/2024.04.18.590027v3)).

## Installing dependencies.

The dependencies are listed in ./requirement.txt.
The key requirements are:
- [looplib](https://github.com/open2c/looplib)
- [polykit](https://github.com/open2c/polykit)
- [hoomd](https://github.com/glotzerlab/hoomd-blue)
- [gsd](https://github.com/glotzerlab/gsd)

## Running

The key script `hoomd_bottlebrush_universal_winding_invariant.py` can be found in the `src/` folder.
To obtain the three published "flagship", run the following commands: 

### SMC3-CAPH-depleted (condensin II-only) mitotic chromosomes , t=30min

```bash
python ./hoomd_bottlebrush_universal_winding_invariant.py cylindrical_helix \
    --rep_A 3.0 \
    --n 500000 \
    --inner_loop_kb 400 \
    --inner_loop_spacing 8 \
    --inner_loops_in_outer 0 \
    --root_loop_spacers \
    --loop_bond_len 1 \
    --nucleosome_box_size 1.64 \
    --bb_angle_wiggle 0.0 \
    --helix_pitch 35 \
    --helix_period_mb 17.0 \
    --fix_terminal_anglepins \
    --fix_turn_anglepins\
    --num_blocks 3000 \
    --replicate 0 \
    --out_folder ../best_models/SMC3-CAPH/ 
```

### SMC3-CAPH2-depleted (condensin I-only) mitotic chromosomes, t=30min

```bash
python ./hoomd_bottlebrush_universal_winding_invariant.py periodic \
    --rep_A 3.0 \
    --n 500000 \
    --inner_loop_kb 100 \
    --inner_loop_spacing 1 \
    --inner_loops_in_outer 0 \
    --root_loop_spacers \
    --loop_bond_len 1 \
    --linear_density 25 \
    --nucleosome_box_size 1.37 \
    --num_blocks 3000 \
    --replicate 0 \
    --out_folder ../best_models/SMC3-CAPH2/ 
```


### SMC3-depleted (condensin I+II) mitotic chromosomes, t=30min
Two set of parameters match the experimental Hi-C equally well:

```bash
python ./hoomd_bottlebrush_universal_winding_invariant.py cylindrical_helix \
    --rep_A 3.0 \
    --n 500000 \
    --inner_loop_kb 100 \
    --inner_loop_spacing 2 \
    --inner_loops_in_outer 4 \
    --outer_loop_spacing 4 \
    --root_loop_spacers \
    --loop_bond_len 1 \
    --nucleosome_box_size 1.37 \
    --bb_angle_wiggle 0.0 \
    --helix_pitch 35 \
    --helix_period_mb 10.0 \
    --fix_terminal_anglepins \
    --fix_turn_anglepins\
    --num_blocks 3000 \
    --replicate 0 \
    --out_folder ../best_models/SMC3_v1/ 
```

or

```bash
python ./hoomd_bottlebrush_universal_winding_invariant.py cylindrical_helix \
    --rep_A 3.0 \
    --n 500000 \
    --inner_loop_kb 100 \
    --inner_loop_spacing 2 \
    --inner_loops_in_outer 4 \
    --outer_loop_spacing 8 \
    --root_loop_spacers \
    --loop_bond_len 1 \
    --nucleosome_box_size 1.37 \
    --bb_angle_wiggle 0.0 \
    --helix_pitch 40 \
    --helix_period_mb 10.0 \
    --fix_terminal_anglepins \
    --fix_turn_anglepins\
    --num_blocks 3000 \
    --replicate 0 \
    --out_folder ../best_models/SMC3_v2/ 
```

## Pre-computed models.

For your convinience, we also provide multiple randomly-generated replicates of each "flagship" model.
These models are stored in the folder `./best_models/` in files named `last_frame.gsd`

