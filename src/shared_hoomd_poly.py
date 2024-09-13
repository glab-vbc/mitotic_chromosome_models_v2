from contextlib import contextmanager
import functools

import numpy as np

import gsd
import gsd.hoomd

import hoomd
import hoomd.md
import hoomd.md.integrate

import looplib
import looplib.random_loop_arrays

import polykit
import polykit.generators.initial_conformations


def _extend_array_or_none(
        array,
        N, 
        default_value,
        dtype,
        second_dim_size=None,
):
    shape = (N, second_dim_size) if second_dim_size else (N,)

    new_chunk = np.full(shape, default_value, dtype)
    if array is None:
        array = new_chunk
    else:
        array = np.concatenate([array, new_chunk])

    return array


def _add_type(snap_property, new_type: str):
    snap_property.types = [] if snap_property.types is None else snap_property.types
    if new_type not in snap_property.types:
        snap_property.types = snap_property.types + [new_type]

    type_id = snap_property.types.index(new_type)
    return type_id


def _extend_typeids(
        snap_property,
        N: int,
        new_type: str,
        ):
    
    # initialize arrays of particle properties if they are missing
    type_id = _add_type(snap_property, new_type)
    
    snap_property.typeid = _extend_array_or_none(snap_property.typeid, N, type_id, np.int32)


def _extend_group(
        snap_property,
        N: int,
):
    snap_property.group = _extend_array_or_none(
        snap_property.group, N, -1, np.int32, second_dim_size=snap_property.M)


def _get_type_id(
    snap_property,
    type_str: str):

    if snap_property.types is None:
        return None
    else:
        return snap_property.types.index(type_str)


def _get_property_idxs_by_id(
    snap_property,
    type_str: str):
    
    return np.where(snap_property.typeid == _get_type_id(snap_property, type_str))[0]


def _add_grouplike_property(
        snap_property,
        group_chunk,
        new_type:str):
    """_summary_

    Args:
        snap_property (_type_): snapshot.bonds, snapshot.dihedrals, etc.
        group_chunk (_type_): _description_
        new_type (str): _description_
    """
    
    snap_property.N += group_chunk.shape[0]
    _extend_typeids(snap_property, group_chunk.shape[0], new_type)
    
    _extend_group(snap_property, group_chunk.shape[0])
    snap_property.group[snap_property.N-group_chunk.shape[0]:] = group_chunk
    
    snap_property.validate()


def add_bonds(
        snap:gsd.hoomd.Frame,
        bonds,
        bond_type:str,
        ):

    _add_grouplike_property(snap.bonds, bonds, bond_type)


def add_angles(
        snap:gsd.hoomd.Frame,
        angles,
        angle_type:str,
        ):

    _add_grouplike_property(snap.angles, angles, angle_type)



def add_dihedrals(
        snap:gsd.hoomd.Frame,
        dihedrals,
        dihedral_type:str,
        ):

    _add_grouplike_property(snap.dihedrals, dihedrals, dihedral_type)


def add_particles(
    snap: gsd.hoomd.Frame,
    N: int,
    particle_type: str,
    positions=None,
    mass=1,
    charge=0,
    diameter=1,
    ):
    
    snap.particles.N += N
    _extend_typeids(snap.particles, N, particle_type)

    snap.particles.position = np.empty(shape=[0, 3], dtype=np.float32) if snap.particles.position is None else snap.particles.position
    if positions is not None:
        snap.particles.position = np.concatenate([snap.particles.position, positions])
    else:
        snap.particles.position = np.concatenate([snap.particles.position, np.zeros(shape=[N, 3], dtype=np.float32)])
    
    snap.particles.mass = _extend_array_or_none(snap.particles.mass, N, mass, np.float32)
    snap.particles.charge = _extend_array_or_none(snap.particles.charge, N, charge, np.float32)
    snap.particles.diameter = _extend_array_or_none(snap.particles.diameter, N, diameter, np.float32)

    snap.particles.validate()




def add_chains(
        snap: gsd.hoomd.Frame,
        chains,
        bond_type: str = 'chain_bond' ,
        angle_type: str = 'chain_angle',
        ):

    for chain in chains:
        chain_bonds = np.vstack([np.arange(chain[0], chain[1]-1), np.arange(chain[0]+1, chain[1])]).T
        add_bonds(snap, chain_bonds, bond_type=bond_type)
        
        if angle_type:
            chain_angles = np.vstack([np.arange(chain[0], chain[1]-2), 
                                      np.arange(chain[0]+1, chain[1]-1), 
                                      np.arange(chain[0]+2, chain[1])]).T
            add_angles(snap, chain_angles, angle_type=angle_type)


def add_exponential_loops(
        snap: gsd.hoomd.Frame,
        chains,
        inner_loop_size,
        inner_loop_spacing,
        loop_spacing_distr='uniform',
        loop_type: str='root_loop',
        ):
    
    for chain in chains:
        n=chain[1]-chain[0]
        loops = looplib.random_loop_arrays.exponential_loop_array(
                N=n, 
                loop_size=inner_loop_size, 
                spacing=inner_loop_spacing, 
                min_loop_size=3,
                loop_spacing_distr=loop_spacing_distr)
        loops += chain[0]

    add_bonds(snap, loops, bond_type=loop_type)


def add_two_layer_exponential_loops(
        snap: gsd.hoomd.Frame,
        chains,
        outer_loop_size,
        outer_loop_spacing,
        inner_loop_size,
        inner_loop_spacing,
        inner_outer_offset,
        outer_loop_type:str='root_loop',
        inner_loop_type:str='nested_loop',
        ):
    
    for chain in chains:
        n=chain[1]-chain[0]
        outer_loops, inner_loops = looplib.random_loop_arrays.two_layer_exponential_loops(
            N=n, 
            outer_loop_size=outer_loop_size, 
            outer_loop_spacing=outer_loop_spacing,
            inner_loop_size=inner_loop_size, 
            inner_loop_spacing=inner_loop_spacing,
            offset=inner_outer_offset)
        inner_loops += chain[0]
        outer_loops += chain[0]

        add_bonds(snap, outer_loops, bond_type=outer_loop_type)
        add_bonds(snap, inner_loops, bond_type=inner_loop_type)    


def add_loop_spacers(
        snap: gsd.hoomd.Frame,
        chains,
        loop_type: str='root_loop',
        spacer_type: str='loop_spacer',
        ):
    """
    TODO: add checks preventing adding spacers between loops on different chains
    """    

    spaced_loops = snap.bonds.group[snap.bonds.typeid == snap.bonds.types.index(loop_type)]    
    spacers = np.vstack([spaced_loops[:-1,1], spaced_loops[1:, 0]]).T

    add_bonds(snap, spacers, bond_type=spacer_type)



def kit_loop_base_angular_pinning(
    snap: gsd.hoomd.Frame,
    forces: dict,

    helix_period_particles: float,
    pin_radius: float=1,
    pin_height: float=1,
    angle_wiggle: float=0.1,
    loop_type_to_pin: str='root_loop',
    dihedral_type: str='anglepin',
    pin_particle_type: str='pin',
    axial_particle_type: str='axis',
    force_name: str='dihedral_pinning',
):
    loops_to_pin = snap.bonds.group[snap.bonds.typeid == snap.bonds.types.index(loop_type_to_pin)]
    n_pins = loops_to_pin.shape[0]
    add_particles(snap, n_pins, pin_particle_type)
    add_particles(snap, 2, axial_particle_type)

    pin_particles_idxs = _get_property_idxs_by_id(snap.particles, pin_particle_type)
    axial_particle_indices = _get_property_idxs_by_id(snap.particles, axial_particle_type)
    pinned_particles_idxs = loops_to_pin[:,0]
        
    new_dihedrals = np.vstack([
        pinned_particles_idxs,
        np.full(n_pins, axial_particle_indices[0]),
        np.full(n_pins, axial_particle_indices[1]),
        pin_particles_idxs
    ]).T

    add_dihedrals(
        snap,
        new_dihedrals,
        dihedral_type=dihedral_type,
        )

    pin_angles = pinned_particles_idxs / helix_period_particles * 2 * np.pi
    snap.particles.position[pin_particles_idxs] = np.vstack([
        np.cos(pin_angles) * pin_radius,
        np.sin(pin_angles) * pin_radius,
        np.full(n_pins, pin_height)]).T

    snap.particles.position[axial_particle_indices[0]] = [0, 0, -pin_height]
    snap.particles.position[axial_particle_indices[1]] = [0, 0, pin_height]

    forces[force_name] = hoomd.md.dihedral.Periodic()
    k = (1/angle_wiggle)**2 if angle_wiggle else 1e-20
    forces[force_name].params[dihedral_type] = dict(k=k, d=-1, n=1, phi0=0)


def add_chain_tip_charge(
        snap: gsd.hoomd.Frame,
        chains,
        ):
    
    for chain in chains:
        snap.particles.charge[chain[0]] = -1
        snap.particles.charge[chain[1]-1] = 1


def kit_chain_stretching_force(
    snap: gsd.hoomd.Frame,
    forces: dict,
    chains,
    strength: float=1.0):
    # this can be re-written to use constant force
    add_chain_tip_charge(snap, chains)

    forces['tip_pinning'] = hoomd.md.external.field.Electric()
    forces['tip_pinning'].E.default = (0,0,strength)



def add_chain_stretching_pins(
    snap: gsd.hoomd.Frame,
    chains,
    linear_density: float,
    pin_particle_type: str='stretching_pin'):
    """
    Stretch chains between their tips by pinning their ends to fixed points in space.
    Changes the type of the first and last particle in each chain to pin_particle_type, 
    which, when used with proper particle filtering, excludes them from simulations.

    Args:
        chains (tuple): tuple of tuples of chain start and end indices
        snap (gsd.hoomd.Frame): a snapshot to modify
        forces (dict): a dictionary of forces to add to
        linear_density (float): linear density of the stretched chains
        pin_particle_type (str, optional): the type of the pin particles. Defaults to 'stretching_pin'.
    """

    for chain in chains:
        n = chain[1] - chain[0]
        length = n / linear_density
        pin_type_id = _add_type(snap.particles, pin_particle_type)
        snap.particles.typeid[chain[0]] = pin_type_id
        snap.particles.typeid[chain[1]-1] = pin_type_id
        snap.particles.position[chain[0]] = [0,0,-length/2]
        snap.particles.position[chain[1]-1] = [0,0,length/2]
        print(f"Chain {chain} stretched to length {length}")


def force_cylindrical_confinement(
    radius: float,
    length: float,
    sigma: float=0.5,
    epsilon: float=1.0,
):
    # U(r) = eps * exp( - 1/2* ((r/sigma)^2 )
    # r(U=1) = sqrt(-2*ln(1/eps)) * sigma
    # At ε=3 and σ=1, r(U=1) = 1.48
    

    energy_at_boundary = 0.3
    wall_shift = np.sqrt(-2*np.log(energy_at_boundary / epsilon)) * sigma

    walls_caps = [
        hoomd.wall.Plane([0,0,-length/2-wall_shift], (0,0,1)),
        hoomd.wall.Plane([0,0, length/2+wall_shift], (0,0,-1)),
    ]

    walls_sides = [
        hoomd.wall.Cylinder(radius=radius+wall_shift, axis=(0,0,1)),
    ]

    # if r_extrap=1, then
    # dU/dr|r=1 = ε * exp(-1/2 * (1/σ)^2) * (-1/σ^2)
    # at sigma=1, dU/dr|r=1 = -ε * exp(-1/2) = -ε * 0.6065


    cyl_wall_force = hoomd.md.external.wall.Gaussian(walls=walls_sides)
    cyl_wall_force.params.default = {"epsilon": epsilon, "sigma": sigma, "r_cut": 3.0*sigma, "r_extrap":1.0}

    cyl_caps_force = hoomd.md.external.wall.Gaussian(walls=walls_caps)
    cyl_caps_force.params.default = {"epsilon": epsilon, "sigma": sigma, "r_cut": 3.0*sigma, "r_extrap":1.0}

    return cyl_wall_force, cyl_caps_force



def force_spherical_confinement(
    radius: float,
    sigma: float=0.5,
    epsilon: float=1.0,
):
    # U(r) = eps * exp( - 1/2* ((r/sigma)^2 )
    # r(U=1) = sqrt(-2*ln(1/eps)) * sigma
    # At ε=3 and σ=1, r(U=1) = 1.48
    
    energy_at_boundary = 0.3
    wall_shift = np.sqrt(-2*np.log(energy_at_boundary / epsilon)) * sigma

    wall = hoomd.wall.Sphere(radius=radius+wall_shift)

    # if r_extrap=1, then
    # dU/dr|r=1 = ε * exp(-1/2 * (1/σ)^2) * (-1/σ^2)
    # at sigma=1, dU/dr|r=1 = -ε * exp(-1/2) = -ε * 0.6065

    sphere_wall_force = hoomd.md.external.wall.Gaussian(walls=[wall])
    sphere_wall_force.params.default = {"epsilon": epsilon, "sigma": sigma, "r_cut": 3.0*sigma, "r_extrap":1.0}

    return sphere_wall_force



def force_cylindrical_axial_exclusion(
    radius: float,
    sigma: float=1.0,
    epsilon: float=2.0,
):
    # U(r) = eps * exp( - 1/2* ((r/sigma)^2 )
    # r(U=1) = sqrt(-2*ln(1/eps)) * sigma
    # At ε=3 and σ=1, r(U=1) = 1.48
    
    walls_sides = [
        hoomd.wall.Cylinder(radius=radius, axis=(0,0,1), inside=False),
    ]

    # if r_extrap=1, then
    # dU/dr|r=1 = ε * exp(-1/2 * (1/σ)^2) * (-1/σ^2)
    # at sigma=1, dU/dr|r=1 = -ε * exp(-1/2) = -ε * 0.6065

    cyl_wall_force = hoomd.md.external.wall.Gaussian(walls=walls_sides)
    cyl_wall_force.params.default = {"epsilon": epsilon, "sigma": sigma, "r_cut": 3.0*sigma, "r_extrap":1.0}

    return cyl_wall_force


def _tag_linker_particles(
    snap: gsd.hoomd.Frame,
    chains,
    loop_type: str='root_loop',
    linker_type: str='linker_type'
    ):

    linker_type_id = _add_type(snap.particles, linker_type)

    loops = snap.bonds.group[snap.bonds.typeid == snap.bonds.types.index(loop_type)]
    for chain in chains:
        chain_loops = loops[(loops[:,0] >= chain[0]) & (loops[:,1] <= chain[1])]
        linkers = np.vstack([chain_loops[:-1,1], chain_loops[1:, 0]]).T
        linkers = np.vstack([[0,chain_loops[0,0]], linkers, [chain_loops[-1,1], chain[-1]-1]])
        for linker in linkers:
            linker_idxs = np.arange(linker[0], linker[1]+1)
            snap.particles.typeid[linker_idxs] = linker_type_id

    return linker_type_id


def tag_terminal_anglepins(
    snap: gsd.hoomd.Frame,
    terminal_dihedral_type: str='terminal_anglepin',
):
    
    type_id = _add_type(snap.dihedrals, terminal_dihedral_type)
    snap.dihedrals.typeid[0] = type_id
    snap.dihedrals.typeid[snap.dihedrals.N-1] = type_id 



def tag_turn_anglepins(
    snap: gsd.hoomd.Frame,
    period_particles: float,
    anglepin_type: str= 'anglepin',
    turn_dihedral_type: str='turn_anglepin',
):

    anglepin_idxs = _get_property_idxs_by_id(snap.dihedrals, anglepin_type)
    pinned_particle_idxs = snap.dihedrals.group[anglepin_idxs][:,0]
    max_pinned_particle_idx = np.max(pinned_particle_idxs)
    ideal_turn_positions = np.arange(0, max_pinned_particle_idx, period_particles)
    real_turn_positions = np.searchsorted(pinned_particle_idxs, ideal_turn_positions)

    type_id = _add_type(snap.dihedrals, turn_dihedral_type)
    snap.dihedrals.typeid[real_turn_positions] = type_id
    

def kit_exclude_linkers_from_axis(
    snap: gsd.hoomd.Frame,
    forces: dict,
    chains,
    radius,
    sigma,
    epsilon,
    loop_type: str='root_loop',
    linker_type: str='linker',
    force_name: str='axial_linker_exclusion',
):
    
    _tag_linker_particles(
        snap,
        chains,
        loop_type=loop_type,
        linker_type=linker_type
        )

    excl_force = force_cylindrical_axial_exclusion(
        radius=radius,
        sigma=1,
        epsilon=0,
    )

    excl_force.params[linker_type] = dict(epsilon=epsilon, sigma=sigma, r_cut=3*sigma, r_extrap=1.0)

    forces[force_name] = excl_force



def init_positions_spacefilling(
        snap: gsd.hoomd.Frame,
        chains: list,
        width: float,
        length: float,
        resize_factor: float=0.8,
        mean_vertical_shift: bool=True,
        ):
    
    init_point = (-width/2, -width/2, -length/2)
    init_positions = []

    for chain in chains:
        chain_init_positions = polykit.generators.initial_conformations.space_filling_curve(
            snap.particles.N, 
            lambda x: (np.abs(x[0]) <= width/2) & (np.abs(x[1]) <= width/2) & (np.abs(x[2]) <= length/2),
            init_point = init_point,
            )
        
        chain_init_positions *= resize_factor
        init_positions.append(chain_init_positions)
        init_point = chain_init_positions[-1]

    if mean_vertical_shift:
        mean_vertical_shift = np.mean(np.concatenate([arr[:,2] for arr in init_positions]))
        for init_positions_chain in init_positions:
            init_positions_chain[:,2] -= mean_vertical_shift
        
    for chain, init_positions_chain in zip(chains, init_positions):
        snap.particles.position[chain[0]:chain[1]] = init_positions_chain
    


def init_positions_rw(
        snap: gsd.hoomd.Frame,
        chains: list,
        width: float,
        length: float,
        radius: float=None,
        resize_factor: float=0.8,
        mean_vertical_shift: bool=True,
    ):
    
    if radius is None and (width is None or length is None):
        raise ValueError("Either radius or width and length must be provided")
    
    init_point = (0, 0, 0)
    init_positions = []

    if radius is None:
        conf_fn = lambda x: (np.abs(x[0]) <= width/2) & (np.abs(x[1]) <= width/2) & (np.abs(x[2]) <= length/2)
    else:
        conf_fn = lambda x: (x[0]**2 + x[1]**2 + x[2]**2 <= radius**2) 

    for chain in chains:
        chain_init_positions = polykit.generators.initial_conformations.create_constrained_random_walk(
            chain[1]-chain[0], 
            conf_fn,
            starting_point = init_point,
            )
        
        chain_init_positions *= resize_factor
        init_positions.append(chain_init_positions)
        init_point = chain_init_positions[-1]

    if mean_vertical_shift:
        mean_vertical_shift = np.mean(np.concatenate([arr[:,2] for arr in init_positions]))
        for init_positions_chain in init_positions:
            init_positions_chain[:,2] -= mean_vertical_shift
        
    for chain, init_positions_chain in zip(chains, init_positions):
        snap.particles.position[chain[0]:chain[1]] = init_positions_chain
    


def init_periodic_conf_polymer(
        snap: gsd.hoomd.Frame,
        chains,
        box: tuple,
        resize_factor: float=0.8,
        mean_vertical_shift: bool=True,
        conformation: str='spacefilling'
    ):

    snap.configuration.box = list(box) + [0, 0, 0]
    if conformation == 'spacefilling':
        init_positions_spacefilling(
                snap,
                chains,
                width=box[0],
                length=box[2],
                resize_factor=resize_factor,
                mean_vertical_shift=mean_vertical_shift,
        )
    elif conformation == 'random_walk':
        init_positions_rw(
                snap,
                chains,
                width=box[0],
                length=box[2],
                resize_factor=resize_factor,
                mean_vertical_shift=mean_vertical_shift,
        )
    else:
        raise NotImplementedError(f"Unsupported conformation type: {conformation}")


def kit_init_cyl_conf_polymer(
        snap: gsd.hoomd.Frame,
        forces: dict,
        chains,
        radius,
        length,
        wall_epsilon,
        ):

    # Periodic boundary conditions by default.
    init_width = 2 * radius / np.sqrt(2)
    snap.configuration.box = [3.0*radius, 3.0*radius, 1.5*length, 0, 0, 0]

    init_positions_spacefilling(
            snap,
            chains,
            width=init_width,
            length=length,
            resize_factor=0.8,
            mean_vertical_shift=True,
    )

    forces['cyl_caps_force'], forces['cyl_wall_force'] = (
        force_cylindrical_confinement(
            radius=radius,
            length=length,
            sigma=0.5,
            epsilon=wall_epsilon,      
        )
    )



def kit_init_spherical_conf_polymer(
        snap: gsd.hoomd.Frame,
        forces: dict,
        chains,
        radius,
        wall_epsilon,
        ):

    # Periodic boundary conditions by default.
    snap.configuration.box = [3.0*radius, 3.0*radius, 3.0*radius, 0, 0, 0]

    init_positions_rw(
            snap,
            chains,
            width=None,
            length=None,
            radius=radius,
            resize_factor=0.8,
            mean_vertical_shift=False,
        )


    forces['sphere_wall_force'] = (
        force_spherical_confinement(
            radius=radius,
            sigma=0.5,
            epsilon=wall_epsilon,      
        )
    )



def force_pin_loopmids_to_side(
    particle_ids_to_pin: list,
    strength: float=5.0,
    direction: tuple=(1,0,0),
    particle_type: str='monomer',
):
    constant = hoomd.md.force.Constant(
        filter=hoomd.filter.Tags(list(particle_ids_to_pin)),

    )
    constant.constant_force[particle_type] = np.array(direction) * strength

    return constant


@contextmanager
def add_gsd_writer(sim, block_size, gsd_path):
    gsd_writer = hoomd.write.GSD(
        filename=gsd_path,
        trigger=hoomd.trigger.Periodic(block_size),
        mode='wb')
    
    sim.operations.writers.append(gsd_writer)
    yield
    sim.operations.writers.pop()


@contextmanager
def add_thermo_tablelog(
        sim, 
        block_size, 
        particle_filter=hoomd.filter.All(),
        log_into_writers=True,
        forces=None
    ):
    

    logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    logger.add(sim, quantities=['timestep', 'tps'])

    thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=particle_filter)
    sim.operations.computes.append(thermodynamic_properties)
    logger.add(
        thermodynamic_properties, 
        quantities=['kinetic_temperature', 'potential_energy'],
        user_name='')
    
    if forces:
        for name, force in forces.items():
            get_energy = functools.partial(getattr, force, 'energy') # had to use this hack to make logging work
            logger[(f'{name}_energy',)] = (get_energy, 'scalar')
            

    table_writer = hoomd.write.Table(
        trigger=hoomd.trigger.Periodic(block_size),
        logger=logger, 
        delimiter="\t",
        pretty=False,
        max_header_len=5)
    table_writer.trigger = hoomd.trigger.Periodic(block_size)

    if log_into_writers:
        for writer in sim.operations.writers:
            writer.logger = logger

    sim.operations.writers.append(table_writer)
    yield
    sim.operations.writers.pop()

    del logger


def run_fire(
        sim,
        forces,
        particle_filter=hoomd.filter.All(),
        n_blocks=10,
        block_size=1000,
        dt=0.05,
        force_tol=5e-2,
        angmom_tol=5e-2,
        energy_tol=5e-2,
        updaters=None,
        flush=True
):

    force_list = list(forces.values()) if isinstance(forces, dict) else forces

    fire = hoomd.md.minimize.FIRE(
        dt=dt,
        force_tol=force_tol,
        angmom_tol=angmom_tol,
        energy_tol=energy_tol,
        # alpha_start=0.999,
        forces=force_list,
        methods=[hoomd.md.methods.ConstantVolume(filter=particle_filter)])
        
    sim.operations.integrator = fire
    sim.run(0)

    for i in range(n_blocks):
        if updaters:
            for updater in updaters:
                updater(sim, forces)
        sim.run(block_size)
        if flush:
            for writer in sim.operations.writers:
                if hasattr(writer, 'flush'):
                    writer.flush()

    for _ in range(len(fire.forces)):
        fire.forces.pop()

    # cleanup
    sim.operations.integrator = None
    del fire

    # restore the velocities lost during FIRE
    sim.state.thermalize_particle_momenta(filter=particle_filter, kT=1.0)


def run_nve(
        sim,
        forces,
        particle_filter=hoomd.filter.All(),
        n_blocks=1000,
        block_size=10000,
        dt=0.05,
        updaters=None,
        flush=True):

    force_list = list(forces.values()) if isinstance(forces, dict) else forces

    nve_filtered = hoomd.md.methods.ConstantVolume(filter=particle_filter)
    nve_filtered_integrator = hoomd.md.Integrator(
        dt=dt, 
        forces=force_list,
        methods=[nve_filtered], 
        )
    sim.operations.integrator = nve_filtered_integrator

    for _ in range(n_blocks):
        if updaters:
            for updater in updaters:
                updater(sim, forces)
                
        sim.run(block_size)

        if flush:
            for writer in sim.operations.writers:
                if hasattr(writer, 'flush'):
                    writer.flush()

    for _ in range(len(nve_filtered_integrator.forces)):
        nve_filtered_integrator.forces.pop()
     
    sim.operations.integrator = None

    del nve_filtered_integrator


def _run_sim(
        sim,
        forces,
        integrator_f=hoomd.md.Integrator,
        method_f=hoomd.md.methods.ConstantVolume,
        particle_filter=hoomd.filter.All(),
        integrator_kwargs={},
        n_blocks=10,
        block_size=1000,
        updaters=None,
        flush=True
):

    force_list = list(forces.values()) if isinstance(forces, dict) else forces

    integrator = integrator_f(
        forces=force_list,
        methods=[method_f(filter=particle_filter)]
        **integrator_kwargs)
        
    sim.operations.integrator = integrator
    sim.run(0)

    for i in range(n_blocks):
        if updaters:
            for updater in updaters:
                updater(sim, forces)
        sim.run(block_size)
        if flush:
            for writer in sim.operations.writers:
                if hasattr(writer, 'flush'):
                    writer.flush()

    for _ in range(len(fire.forces)):
        integrator.forces.pop()

    # cleanup
    sim.operations.integrator = None
    del fire

    # restore the velocities lost during FIRE
    sim.state.thermalize_particle_momenta(filter=particle_filter, kT=1.0)

