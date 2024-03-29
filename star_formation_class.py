"""
Grouped star formation class. This class is taken from
star_forming_region_class.py in Ekster, which the development will mainly
be in: https://github.com/rieder/ekster

References:
-   Liow K. Y., Rieder S., Dobbs C. L., Jaffa S. E., 2022, MNRAS, 510, 2657
-   Rieder S., Dobbs C. L., Bending T., Liow K. Y., Wurster J., 2022, MNRAS,
    509, 6155
"""

import numpy

from amuse.datamodel import Particles
from amuse.units import units
from amuse.ic.brokenimf import new_kroupa_mass_distribution
from amuse.units.trigo import sin, cos
from amuse.ext.masc.cluster import new_masses


def generate_next_mass(
        lower_mass_limit=0.5|units.MSun,
        upper_mass_limit=100|units.MSun,
):
    "Generate list of masses of next star/stars to form"
    return new_kroupa_mass_distribution(
        1,
        mass_min=lower_mass_limit,
        mass_max=upper_mass_limit,
    )


def assign_sink_group(
        sink,
        sink_particles,
        group_radius=1|units.pc,
        group_speed=0.2|units.kms,
        group_age=0.1|units.Myr,
):
    """
    Assign group index to sink particle. All initialised sinks must
    have a group index.
    """

    if not hasattr(sink, "in_group"):
        sink.in_group = 0

    number_of_groups = sink_particles.in_group.max()
    initialised = sink.initialised or False
    if not initialised:

        # Check if this sink belongs to any existing groups. Must
        # pass all checks.
        smallest_Etot = numpy.inf | units.J
        fail1 = fail2 = fail3 = fail4 = 0
        for i in range(number_of_groups):
            i += 1   # Change to one-based index
            group_i = sink_particles[sink_particles.in_group == i]

            # Check 1: see if this sink is within the sampling radius
            # from the center of mass of the i-th group.
            distance_from_group_com = (
                sink.position - group_i.center_of_mass()
            ).length()
            if distance_from_group_com > group_radius:
                fail1 += 1
                continue

            # Check 2: see if this sink is within the sampling
            # velocity from the center-of-mass velocity of the group
            speed_from_group_com = (
                sink.velocity - group_i.center_of_mass_velocity()
            ).length()
            if speed_from_group_com > group_speed:
                fail2 += 1
                continue

            # Check 3: see if 'the sink' is similar in age with the group
            age_difference = sink.birth_time - group_i.birth_time.min()
            if age_difference > group_age:
                fail3 += 1
                continue

            group_and_sink = group_i.copy()
            group_and_sink.add_particle(sink.copy())
            Etot = (
                group_and_sink.kinetic_energy()
                + group_and_sink.potential_energy()
            )
            # Check 4: see if this sink is the most bound to this
            # group
            if Etot > smallest_Etot:
                fail4 += 1
                continue

            # At this point, this sink passes all checks
            smallest_Etot = Etot
            sink.in_group = i

        # If this sink is still unassigned to any of the groups,
        # create its own group
        if sink.in_group == 0:
            sink.in_group = number_of_groups + 1

        sink.initialised = True

    return sink


def form_stars_from_group(
    group_index,
    sink_particles,
    lower_mass_limit=0.5|units.MSun,
    upper_mass_limit=100|units.MSun,
    local_sound_speed=0.3 | units.kms,
    minimum_sink_mass=0.01 | units.MSun,
    randomseed=None,
    shrink_sinks=True,
    **keyword_arguments
):
    """
    Form stars from specific group of sinks.
    """
    if randomseed is not None:
        numpy.random.seed(randomseed)

    # Consider only group with input group index from here onwards.
    group = sink_particles[sink_particles.in_group == group_index]

    # Sanity check: group must have at least a sink
    if group.is_empty():
        print(
            "ERROR: There is no sink in the group! Check group assignment."
        )
        exit()

    number_of_sinks = len(group)
    group_mass = group.total_mass()

    next_mass = generate_next_mass(
        lower_mass_limit=lower_mass_limit,
        upper_mass_limit=upper_mass_limit,
    )[0][0]
    try:
        # Within a group, group_next_primary_mass values are either
        # a mass, or 0 MSun. If all values are 0 MSun, this is a
        # new group. Else, only interested on the non-zero value. The
        # non-zero values are the same.
        if group.group_next_primary_mass.max() == 0 | units.MSun:
            group.group_next_primary_mass = next_mass
        else:
            next_mass = group.group_next_primary_mass.max()
    # This happens for the first ever assignment of this attribute
    except AttributeError:
        group.group_next_primary_mass = next_mass

    #logger.info("Next mass is %s", next_mass)

    if group_mass < next_mass:
        # Group is not massive enough for next star
        return None

    # Form stars from the leftover group sink mass
    mass_left = group_mass - next_mass
    masses = new_masses(
        stellar_mass=mass_left,
        lower_mass_limit=lower_mass_limit,
        upper_mass_limit=upper_mass_limit,
        initial_mass_function="kroupa",
        exceed_mass=False,
    )
    number_of_stars = len(masses)

    if number_of_stars == 0:
        # No stars created for this group
        return None

    new_stars = Particles(number_of_stars)
    new_stars.age = 0 | units.Myr
    new_stars[0].mass = next_mass
    new_stars[1:].mass = masses[:-1]
    group.group_next_primary_mass = masses[-1]
    new_stars = new_stars.sorted_by_attribute("mass").reversed()
    new_stars.in_group = group_index

    # Create placeholders for attributes of new_stars
    new_stars.position = [0, 0, 0] | units.pc
    new_stars.velocity = [0, 0, 0] | units.kms
    new_stars.origin_cloud = group[0].key
    new_stars.star_forming_radius = 0 | units.pc
    new_stars.star_forming_u = local_sound_speed**2

    # Don't mess with the actual group sink particle set.
    star_forming_regions = group.copy()
    star_forming_regions.sorted_by_attribute("mass").reversed()

    # Generate a probability list of star forming region indices the
    # stars should associate to
    probabilities = (
        star_forming_regions.mass/star_forming_regions.mass.sum()
    )
    probabilities /= probabilities.sum()    # Ensure sum is exactly 1

    # Create index list of star forming regions from probability list
    sample = numpy.random.choice(
        len(star_forming_regions), number_of_stars, p=probabilities
    )

    # Assign the stars to the sampled star forming regions
    star_forming_regions_sampled = star_forming_regions[sample]
    new_stars.position = star_forming_regions_sampled.position
    new_stars.velocity = star_forming_regions_sampled.velocity
    new_stars.origin_cloud = star_forming_regions_sampled.key
    new_stars.star_forming_radius = star_forming_regions_sampled.radius
    try:
        new_stars.star_forming_u = star_forming_regions_sampled.u
    except AttributeError:
        new_stars.star_forming_u = local_sound_speed**2

    # Random position of stars within the sink radius they assigned to
    rho = (
        numpy.random.random(number_of_stars) * new_stars.star_forming_radius
    )
    theta = (
        numpy.random.random(number_of_stars)
        * (2 * numpy.pi | units.rad)
    )
    phi = (
        numpy.random.random(number_of_stars) * numpy.pi | units.rad
    )
    x = (rho * sin(phi) * cos(theta)).value_in(units.pc)
    y = (rho * sin(phi) * sin(theta)).value_in(units.pc)
    z = (rho * cos(phi)).value_in(units.pc)

    dX = list(zip(*[x, y, z])) | units.pc

    # Random velocity, sample magnitude from gaussian with local sound speed
    # like Wall et al (2019)
    # temperature = 10 | units.K

    # or (gamma * local_pressure / density).sqrt()
    velocity_magnitude = numpy.random.normal(
        # loc=0.0,  # <- since we already added the velocity of the sink
        scale=new_stars.star_forming_u.sqrt().value_in(units.kms),
        size=number_of_stars,
    ) | units.kms
    velocity_theta = (
        numpy.random.random(number_of_stars)
        * (2 * numpy.pi | units.rad)
    )
    velocity_phi = (
        numpy.random.random(number_of_stars)
        * (numpy.pi | units.rad)
    )
    vx = (
        velocity_magnitude * sin(velocity_phi) * cos(velocity_theta)
    ).value_in(units.kms)
    vy = (
        velocity_magnitude * sin(velocity_phi) * sin(velocity_theta)
    ).value_in(units.kms)
    vz = (
        velocity_magnitude * cos(velocity_phi)
    ).value_in(units.kms)

    dV = list(zip(*[vx, vy, vz])) | units.kms

    #logger.info("Updating new stars...")
    new_stars.position += dX
    new_stars.velocity += dV

    # For Pentacle, this is the PP radius
    new_stars.radius = 0.05 | units.parsec

    # Remove sink mass according to the position of stars
    excess_star_mass = 0 | units.MSun
    for s in group:
        #logger.info('Sink mass before reduction: %s', s.mass.in_(units.MSun))
        total_star_mass_nearby = (
            new_stars[new_stars.origin_cloud == s.key]
        ).total_mass()

        # To prevent sink mass becomes negative
        if s.mass > minimum_sink_mass:
            if (s.mass - total_star_mass_nearby) <= minimum_sink_mass:
                excess_star_mass += (
                    total_star_mass_nearby - s.mass + minimum_sink_mass
                )
                #logger.info(
                #    'Sink mass goes below %s; excess mass is now %s',
                #    minimum_sink_mass.in_(units.MSun),
                #    excess_star_mass.in_(units.MSun)
                #)
                s.mass = minimum_sink_mass
            else:
                s.mass -= total_star_mass_nearby
        else:
            excess_star_mass += total_star_mass_nearby
            #logger.info(
            #    'Sink mass is already <= minimum mass allowed; '
            #    'excess mass is now %s',
            #    excess_star_mass.in_(units.MSun)
            #)

        #logger.info('Sink mass after reduction: %s', s.mass.in_(units.MSun))

    # Reduce all sinks in group equally with the excess star mass
    #logger.info('Reducing all sink mass equally with excess star mass...')
    mass_ratio = 1 - excess_star_mass/group.total_mass()
    group.mass *= mass_ratio
    #
    # logger.info(
    #     "Total sink mass in group after sink mass reduction: %s",
    #     group.total_mass().in_(units.MSun)
    # )

    if shrink_sinks:
        group.radius = (
            (group.mass / group.initial_density)
            / (4/3 * numpy.pi)
        )**(1/3)

    return new_stars
