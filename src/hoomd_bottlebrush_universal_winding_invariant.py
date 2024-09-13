import argparse
import pathlib
import pprint
import logging
import functools
import types

import numpy as np

import hoomd
import hoomd.md
import hoomd.md.integrate
import gsd
import gsd.hoomd

import shared_hoomd_poly

print = functools.partial(print, flush=True)

parser = argparse.ArgumentParser()

subparsers = parser.add_subparsers(help='Style of confinement', dest='conf_type')


parser_periodic = subparsers.add_parser('periodic', help='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_cyl = subparsers.add_parser('cylindrical', help='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_helix = subparsers.add_parser('cylindrical_helix', help='',formatter_class=argparse.ArgumentDefaultsHelpFormatter)


all_parsers = [parser_cyl, parser_helix, parser_periodic]

for cur_parser in all_parsers:   
    size_arg_group = cur_parser.add_mutually_exclusive_group(required=True)
    size_arg_group.add_argument(
        '--n', 
        type=int, 
        default=0, 
        help="The number of particles in the system. Provide either --n or --loop_n.")
    size_arg_group.add_argument(
        '--loop_n', 
        type=int, 
        default=0, 
        help="The number of loops in the system. Provide either --n or --loop_n.")

    bb_arg_group = cur_parser.add_argument_group(title='Bottlebrush', description=None)
    bb_arg_group.add_argument(
        '--inner_loop_kb', 
        type=float, 
        default=100,
        help="The average size of loops in kb.")
    bb_arg_group.add_argument(
        '--inner_loops_in_outer', 
        type=float, 
        default=0,
        help="The average size of loops in kb.")
    bb_arg_group.add_argument(
        '--inner_loop_spacing', 
        type=int, 
        default=5, 
        help="The size of a linker between two consecutive loops.")
    bb_arg_group.add_argument(
        '--outer_loop_spacing', 
        type=int, 
        default=6, 
        help="The size of a linker between two consecutive loops.")
    bb_arg_group.add_argument(
        '--inner_outer_offset', 
        type=int, 
        default=5, 
        help="The size of a linker between two consecutive loops.")
    bb_arg_group.add_argument(
        '--loop_bond_len',
        type=float, 
        default=1.0, 
        help="The length of bonds supporting loops.")
    bb_arg_group.add_argument(
        '--root_loop_spacers',
        action='store_true',
        help='If provided, the spacers between root loops will be straight.')

    conf_arg_group = cur_parser.add_argument_group(title='Confinement', description=None)
    conf_arg_group.add_argument(
        '--nucleosome_box_size', 
        type=float, 
        default=1.4, 
        help="The specific volume of chromosome per nucleosome, provided as a cubic root.")

    chain_force_arg_group = cur_parser.add_argument_group(title='Chain forces', description=None)
    chain_force_arg_group.add_argument(
        '--rep_A',
        type=float, 
        default=3.0, 
        help="The maximal energy of repulsion between particles.")
    chain_force_arg_group.add_argument(
        '--bond_k', 
        type=float, 
        default=10.0, 
        help="The stiffness of all harmonic bonds (polymer, loops, spacers).")

    integration_arg_group = cur_parser.add_argument_group(title='Integration', description=None)
    integration_arg_group.add_argument(
        '--num_blocks', 
        type=int, 
        default=1000, 
        help="The number of blocks to simulate.")
    integration_arg_group.add_argument(
        '--dt',
        type=float, 
        default=0.05, 
        help="The timestep of integration.")
    integration_arg_group.add_argument(
        '--gamma',
        type=float, 
        default=1.0, 
        help="The DPD friction coefficient.")


    misc_arg_group = cur_parser.add_argument_group(title='Miscellaneous', description=None)
    misc_arg_group.add_argument(
        '--replicate', 
        type=int, 
        default=0, 
        help="The replicate index.")
    misc_arg_group.add_argument(
        '--init_conf', 
        type=str, 
        choices=['spacefilling', 'random_walk'],
        default='spacefilling', 
        help="The initial conformation of the polymer. Possible values: 'spacefilling', 'random_walk'.")
    misc_arg_group.add_argument(
        '--out_folder', 
        type=str, 
        default='./', 
        help="The root folder where the results will be stored.")
    misc_arg_group.add_argument(
        '--overwrite', 
        action='store_true')

# Cylinder specific arguments
cyl_arg_group = parser_cyl.add_argument_group(title='Cylindrical confinement', description=None)
cyl_arg_group.add_argument(
    '--linear_density', 
    type=float, 
    default=None, 
    help="The linear density of chromatin along the cylinder in Mb/mum. ")
cyl_arg_group.add_argument(
    '--cyl_conf_epsilon', 
    type=float, 
    default=1.0, 
    help="The wall repulsion strength. If set to 0, no cylindrical confinement is applied.")


helix_arg_group = parser_helix.add_argument_group(title='Helical pinning in cylindrical confinement', description=None)
helix_arg_group.add_argument(
    '--linear_density', 
    type=float, 
    default=None, 
    help="The linear density of chromatin along the cylinder in Mb/mum. Specify only two arguments out of linear_density, helix_period_mb, helix_pitch.")
helix_arg_group.add_argument(
    '--helix_period_mb', 
    type=float, 
    default=10.0, 
    help="Period of the spiral in MB.")
helix_arg_group.add_argument(
    '--helix_pitch', 
    type=float, 
    default=40, 
    help="The pitch of the helix.")

helix_arg_group.add_argument(
    '--bb_angle_wiggle', 
    type=float, 
    default=0.8, 
    help="Angular wiggle of the backbone particles in radian.")
helix_arg_group.add_argument(
    '--fix_terminal_anglepins', 
    action='store_true',
    help="If true, fix the angles of the first and last loop base.")
helix_arg_group.add_argument(
    '--fix_turn_anglepins', 
    action='store_true',
    help="If true, fix the angles of one particle in each turn of the helix.")
helix_arg_group.add_argument(
    '--cyl_conf_epsilon', 
    type=float, 
    default=1.0, 
    help="The wall repulsion strength. If set to 0, no cylindrical confinement is applied.")


# Periodic-box specific arguments
linear_stretching = parser_periodic.add_argument_group(title='Helical pinning in cylindrical confinement', description=None)
linear_stretching.add_argument(
    '--linear_density', 
    type=float, 
    default=None, 
    help="If prodived, the tips of the chain will be pinned to stretch it to the specified linear density in MB/mum.")


args = parser.parse_args()

p = types.SimpleNamespace()
p.BP_PP = 200

p.INNER_LOOP_SIZE_KB = args.inner_loop_kb
p.INNER_LOOP_SIZE = int(np.round(args.inner_loop_kb * 1000.0 / p.BP_PP))
p.INNER_LOOP_SPACING = args.inner_loop_spacing

p.INNER_LOOPS_IN_OUTER = args.inner_loops_in_outer
p.OUTER_LOOP_SIZE = int(np.round(args.inner_loop_kb * 1000.0 / p.BP_PP * p.INNER_LOOPS_IN_OUTER))
p.OUTER_LOOP_SPACING = args.outer_loop_spacing
p.INNER_OUTER_OFFSET = args.inner_outer_offset

p.LOOP_BOND_LEN = args.loop_bond_len
p.ROOT_SPACERS = args.root_loop_spacers

if args.n and not args.loop_n:
    p.N = args.n
    p.LOOP_N = p.N // p.INNER_LOOP_SIZE if p.INNER_LOOP_SIZE else 0
    p.N_CALC_METHOD = 'n'
elif not args.n and args.loop_n:
    p.LOOP_N = args.loop_n
    p.N = p.INNER_LOOP_SIZE*p.LOOP_N
    p.N_CALC_METHOD = 'loop_n'
else:
    raise ValueError('Provide either --n or --loop_n')

p.REP_A = args.rep_A
p.DT = args.dt
p.GAMMA = args.gamma
p.BOND_K = args.bond_k

p.REPLICATE = args.replicate
p.NUM_BLOCKS = args.num_blocks
p.BLOCK_SIZE = 10000
p.FIRE_BLOCK_SIZE = 1000
p.OUT_ROOT_FOLDER = pathlib.Path(args.out_folder)


p.BOX_SIZE_PP = args.nucleosome_box_size
p.V = p.N * (p.BOX_SIZE_PP ** 3)

if args.conf_type == 'cylindrical':
    p.CONFINEMENT_TYPE = 'cylindrical'
    p.LINEAR_DENSITY_MB_MUM = args.linear_density

    p.CYL_CONF_EPSILON = args.cyl_conf_epsilon
    p.LENGTH = p.BP_PP * p.N / 1e6 / p.LINEAR_DENSITY_MB_MUM * 100
    p.RADIUS = (p.V / p.LENGTH / np.pi) ** (1/2)

elif args.conf_type == 'cylindrical_helix':
    p.CONFINEMENT_TYPE = 'cylindrical_helix'
    if sum(not v for v in [args.helix_period_mb, args.helix_pitch, args.linear_density]) != 1:
        raise ValueError('Specify exactly two arguments out of helix_period_mb, helix_pitch, and linear_density.')
    if not args.linear_density:
        p.HELIX_PERIOD_MB = args.helix_period_mb 
        p.HELIX_PITCH = args.helix_pitch
        p.LINEAR_DENSITY_MB_MUM = p.HELIX_PERIOD_MB / (p.HELIX_PITCH / 100)
    elif not args.helix_period_mb:
        p.LINEAR_DENSITY_MB_MUM = args.linear_density
        p.HELIX_PITCH = args.helix_pitch
        p.HELIX_PERIOD_MB = p.LINEAR_DENSITY_MB_MUM * p.HELIX_PITCH / 100
    elif not args.helix_pitch:
        p.LINEAR_DENSITY_MB_MUM = args.linear_density
        p.HELIX_PERIOD_MB = args.helix_period_mb
        p.HELIX_PITCH = p.HELIX_PERIOD_MB / (p.LINEAR_DENSITY_MB_MUM / 100)


    p.HELIX_PERIOD_PARTICLES = p.HELIX_PERIOD_MB * 1e6 / p.BP_PP
    p.LINEAR_DENSITY = p.LINEAR_DENSITY_MB_MUM * 1e6 / p.BP_PP / 100
    p.BB_ANGLE_WIGGLE = args.bb_angle_wiggle
    p.FIX_TERMINAL_ANGLEPINS = args.fix_terminal_anglepins
    p.FIX_TURN_ANGLEPINS = args.fix_turn_anglepins
    
    p.CYL_CONF_EPSILON = args.cyl_conf_epsilon
    p.LENGTH = p.BP_PP * p.N / 1e6 / p.LINEAR_DENSITY_MB_MUM * 100
    p.RADIUS = (p.V / p.LENGTH / np.pi) ** (1/2)
    

elif args.conf_type == 'periodic':
    p.CONFINEMENT_TYPE = 'periodic'
    p.LINEAR_DENSITY_MB_MUM = args.linear_density
    if p.LINEAR_DENSITY_MB_MUM is None:
        p.LENGTH = (p.V) ** (1/3)
        p.LINEAR_DENSITY = None
    else:
        p.LENGTH = p.BP_PP * p.N / 1e6 / p.LINEAR_DENSITY_MB_MUM * 100
        p.LINEAR_DENSITY = p.LINEAR_DENSITY_MB_MUM * 1e6 / p.BP_PP / 100
    p.WIDTH = (p.V / p.LENGTH) ** (1/2)
    p.INIT_CONF = args.init_conf

def make_sim_name(p):
    sim_name_dict = {}
    if p.N_CALC_METHOD == 'loop_n':
        sim_name_dict['LoopN'] = args.loop_n
    elif p.N_CALC_METHOD == 'n':
        sim_name_dict['N'] = args.n
    
    sim_name_dict.update(dict(
        InLoop=p.INNER_LOOP_SIZE_KB,
        OutLoop=p.INNER_LOOPS_IN_OUTER,
        InSpacing=p.INNER_LOOP_SPACING,
        OutSpacing=p.OUTER_LOOP_SPACING,
        IOOffset=p.INNER_OUTER_OFFSET,
        LoopBondLen=p.LOOP_BOND_LEN,
    ))

    if p.ROOT_SPACERS:
        sim_name_dict['RootSpacers'] = 1
    

    if p.CONFINEMENT_TYPE == 'cylindrical':
        sim_name_dict.update(dict(
            ConfType='Cyl',
            LinDens=p.LINEAR_DENSITY_MB_MUM,
            WallE=p.CYL_CONF_EPSILON,
        ))

    elif p.CONFINEMENT_TYPE == 'cylindrical_helix':
        sim_name_dict.update(dict(
            ConfType='CylHel',
            PeriodMb=p.HELIX_PERIOD_MB,
            Pitch=p.HELIX_PITCH,
            BBAngWiggle=p.BB_ANGLE_WIGGLE,
            WallE=p.CYL_CONF_EPSILON,
        ))
        if p.FIX_TERMINAL_ANGLEPINS:
            sim_name_dict['FixTerms'] = 1
        if p.FIX_TURN_ANGLEPINS:
            sim_name_dict['FixTurns'] = 1

    elif p.CONFINEMENT_TYPE == 'periodic':
        sim_name_dict.update(dict(
            ConfType='Periodic',
            LinDens=p.LINEAR_DENSITY_MB_MUM))


    sim_name_dict.update(dict(
        BoxPP=p.BOX_SIZE_PP,
        RepA=p.REP_A,
        G=p.GAMMA,
        DT=p.DT,
        R=p.REPLICATE,
    ))

    return sim_name_dict


    
p.SIM_NAME = "-".join(f'{n}_{v}' for n, v in make_sim_name(p).items())

p.OUT_FOLDER = p.OUT_ROOT_FOLDER / p.SIM_NAME

if p.OUT_FOLDER.exists():
    if not args.overwrite:
        raise Exception('The output folder already exists!')
else:
    p.OUT_FOLDER.mkdir(parents=True)

p.OUT_INIT_PATH = p.OUT_FOLDER/'init.gsd'
p.OUT_OPTIMIZED_PATH = p.OUT_FOLDER/'optimized.gsd'
p.OUT_TRAJ_PATH = p.OUT_FOLDER/'traj.gsd'
p.LOG_PATH = p.OUT_FOLDER/'out.log'


logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(p.LOG_PATH),
        logging.StreamHandler()
    ]
)

logging.debug(pprint.pformat(p.__dict__))


chains = [(0,p.N)]

snapshot = gsd.hoomd.Frame()
forces = {}

shared_hoomd_poly.add_particles(
    snapshot,
    p.N,
    particle_type='monomer',
    )


shared_hoomd_poly.add_chains(
        snapshot,
        chains,
        bond_type='chain_bond',
        )

if p.OUTER_LOOP_SIZE:
    shared_hoomd_poly.add_two_layer_exponential_loops(
        snapshot,
        chains,
        p.OUTER_LOOP_SIZE,
        p.OUTER_LOOP_SPACING,
        p.INNER_LOOP_SIZE,
        p.INNER_LOOP_SPACING,
        p.INNER_OUTER_OFFSET,
        )

elif p.INNER_LOOP_SIZE:
    shared_hoomd_poly.add_exponential_loops(
        snapshot,
        chains,
        p.INNER_LOOP_SIZE,
        p.INNER_LOOP_SPACING,
        )

if p.ROOT_SPACERS:
    shared_hoomd_poly.add_loop_spacers(snapshot, chains)



# creates forces, first bonds
forces['harmonic'] = hoomd.md.bond.Harmonic()

# bonds must be configured for each bond type
forces['harmonic'].params.default = dict(k=p.BOND_K, r0=1.0)

if p.ROOT_SPACERS:
    if p.OUTER_LOOP_SIZE:
        forces['harmonic'].params['loop_spacer'] = dict(k=p.BOND_K, r0=float(p.OUTER_LOOP_SPACING))
    else:
        forces['harmonic'].params['loop_spacer'] = dict(k=p.BOND_K, r0=float(p.INNER_LOOP_SPACING))

if p.CONFINEMENT_TYPE == 'cylindrical' or p.CONFINEMENT_TYPE == 'cylindrical_helix':
    shared_hoomd_poly.kit_init_cyl_conf_polymer(
            snap=snapshot,
            forces=forces,
            chains=chains, 
            radius=p.RADIUS,
            length=p.LENGTH, 
            wall_epsilon=p.CYL_CONF_EPSILON,
            )
elif p.CONFINEMENT_TYPE == 'periodic':
    shared_hoomd_poly.init_periodic_conf_polymer(
            snapshot,
            chains, 
            [p.WIDTH, p.WIDTH, p.LENGTH],
            resize_factor=0.5,
            conformation=p.INIT_CONF,
            )



# nlist = hoomd.md.nlist.Cell(buffer=0.4, exclusions=tuple())
nlist = hoomd.md.nlist.Tree(buffer=0.4)#, exclusions=tuple())
# nlist = hoomd.md.nlist.Stencil(cell_width=10.0, buffer=0.4)

forces['dpd'] = hoomd.md.pair.DPD(nlist=nlist, kT=1.0, default_r_cut=1.0)
forces['dpd'].params.default = dict(A=p.REP_A, gamma=p.GAMMA)

if p.CONFINEMENT_TYPE == 'cylindrical':
    shared_hoomd_poly.kit_chain_stretching_force(snapshot, forces, chains)    

elif p.CONFINEMENT_TYPE == 'cylindrical_helix':
    shared_hoomd_poly.kit_loop_base_angular_pinning(
        snap=snapshot,
        forces=forces,
        helix_period_particles=p.HELIX_PERIOD_PARTICLES,
        angle_wiggle=0.1,
    )

    shared_hoomd_poly.kit_exclude_linkers_from_axis(
        snapshot,
        forces,
        chains,
        radius=2, # this used to be 5, but proved too wide for small inter-loop spacings
        sigma=1,
        epsilon=0)
    
    if p.FIX_TURN_ANGLEPINS:
        shared_hoomd_poly.tag_turn_anglepins(snapshot, p.HELIX_PERIOD_PARTICLES)
        forces['dihedral_pinning'].params['turn_anglepin'] = forces['dihedral_pinning'].params['anglepin']
    if p.FIX_TERMINAL_ANGLEPINS:
        shared_hoomd_poly.tag_terminal_anglepins(snapshot)
        forces['dihedral_pinning'].params['terminal_anglepin'] = forces['dihedral_pinning'].params['anglepin']

    shared_hoomd_poly.kit_chain_stretching_force(snapshot, forces, chains)

elif p.CONFINEMENT_TYPE == 'periodic':
    if p.LINEAR_DENSITY_MB_MUM is not None:
        shared_hoomd_poly.add_chain_stretching_pins(snapshot, chains, p.LINEAR_DENSITY)






with gsd.hoomd.open(name=p.OUT_INIT_PATH, mode='w') as f:
    f.append(snapshot)

gpu = hoomd.device.GPU()
sim = hoomd.Simulation(device=gpu, seed=1)

# set coordinates
sim.create_state_from_gsd(filename=p.OUT_INIT_PATH)

# set speeds
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.0)

if 'linker' in snapshot.particles.types:
    particle_filter = hoomd.filter.Type(['monomer', 'linker'])
else:
    particle_filter = hoomd.filter.Type(['monomer'])

def enable_axial_linker_exclusion(sim, forces, epsilon=6, sigma=1):
    if 'axial_linker_exclusion' in forces:
        if sim.timestep // p.FIRE_BLOCK_SIZE == 10:
            forces['axial_linker_exclusion'].params['linker'] = dict(epsilon=epsilon, sigma=sigma, r_cut=3*sigma, r_extrap=1.0)
            logging.debug(f"Turned axial exclusion on at timestep {sim.timestep}")

fire_updaters = [enable_axial_linker_exclusion] if p.CONFINEMENT_TYPE == 'cylindrical_helix' else []

with shared_hoomd_poly.add_gsd_writer(sim, p.FIRE_BLOCK_SIZE, p.OUT_OPTIMIZED_PATH):
    with shared_hoomd_poly.add_thermo_tablelog(
        sim, p.FIRE_BLOCK_SIZE, particle_filter=particle_filter, log_into_writers=True, forces=forces):

        shared_hoomd_poly.run_fire(
            sim,
            forces,
            particle_filter=particle_filter,
            n_blocks=20,
            block_size=p.FIRE_BLOCK_SIZE,
            updaters=fire_updaters,
            )
        
def update_angular_pinning(sim, forces):
    if 'diheral_pinning' in forces:
        if sim.timestep // p.BLOCK_SIZE == 20:
            new_k = (1/p.BB_ANGLE_WIGGLE)**2 if p.BB_ANGLE_WIGGLE else 1e-20
            forces['dihedral_pinning'].params['anglepin'] = dict(k=new_k, d=-1, n=1, phi0=0)
            logging.debug(f"Update angular pinning on at timestep {sim.timestep}, new_k={new_k}")

main_run_updaters = [update_angular_pinning] if p.CONFINEMENT_TYPE == 'cylindrical_helix' else []

with shared_hoomd_poly.add_gsd_writer(sim, p.BLOCK_SIZE, p.OUT_TRAJ_PATH):
    with shared_hoomd_poly.add_thermo_tablelog(
            sim, p.BLOCK_SIZE, particle_filter=particle_filter, log_into_writers=True, forces=forces):
        
        shared_hoomd_poly.run_nve(
                sim,
                forces,
                particle_filter=particle_filter,
                n_blocks=p.NUM_BLOCKS,
                block_size=p.BLOCK_SIZE,
                dt=p.DT,
                updaters=main_run_updaters,
        )


