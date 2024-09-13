# mitotic_chromosome_models_v2
The scripts and an accompanying library generating polymer models of mitotic chromosomes.   


## SMC3-CAPH-depleted mitotic chromosomes, t=30min

bash```
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


## SMC3-CAPH2-depleted mitotic chromosomes, t=30min

bash```
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


## SMC3-depleted mitotic chromosomes, t=30min

bash```
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

bash```
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