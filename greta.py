import argparse

import numpy

from amuse.io import read_set_from_file, write_set_to_file
from amuse.units import units
from amuse.datamodel import Particles

from star_formation_class import assign_sink_group, form_stars_from_group

def new_argument_parser():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        dest='sinkfilename',
        default='',
        help='AMUSE sink particle set',
    )
    parser.add_argument(
        '-o',
        dest='starfilename',
        default='stars.amuse',
        help='Output name for AMUSE star particle set [stars.amuse]',
    )
    parser.add_argument(
        '-d',
        dest='grouping_distance',
        type=float,
        default=0,
        help='Grouping distance parameter in pc [0]',
    )
    parser.add_argument(
        '-v',
        dest='grouping_speed',
        type=float,
        default=0,
        help='Grouping speed parameter in km/s [0]',
    )
    parser.add_argument(
        '-t',
        dest='grouping_age',
        type=float,
        default=0,
        help='Grouping age parameter in Myr [0]'
    )
    parser.add_argument(
        '--lower-limit',
        dest='mass_lower',
        type=float,
        default=0.5,
        help='Lower star mass limit in MSun [0.5]'
    )
    parser.add_argument(
        '--upper-limit',
        dest='mass_upper',
        type=float,
        default=100,
        help='Upper star mass limit in MSun [100]'
    )
    return parser.parse_args()


def main():
    o = new_argument_parser()
    sinkfilename = o.sinkfilename
    starfilename = o.starfilename
    grouping_distance = o.grouping_distance | units.pc
    grouping_speed = o.grouping_speed | units.kms
    grouping_age = o.grouping_age | units.Myr
    mass_lower = o.mass_lower | units.MSun
    mass_upper = o.mass_upper | units.MSun

    sinks = read_set_from_file(sinkfilename, 'amuse')
    print(f"{len(sinks)} sink particles")

    if not hasattr(sinks, 'position'):
        msg = (
            "ERROR: Sink particles have no positions! Exiting program."
        )
        print(msg)
        exit()

    if not hassattr(sinks, 'velocity'):
        msg = (
            "WARNING: Sink particles have no velocities! Assuming they are "
            + "stationary."
        )
        print(msg)
        sinks.velocity = 0 | units.kms


    if not hasattr(sinks, 'radius'):
        msg = (
            "WARNING: Sink particles have no attribute 'radius', setting to "
            + "0.1 pc. This does not matter because this program is not "
            + "supposed to be used to generate initial conditions."
        )
        print(msg)
        sinks.radius = 0.1 | units.pc

    if not hassattr(sinks, 'birth_time'):
        msg = (
            "Sink particles have no attribute 'birth_time', setting to "
            + "0 Myr, i.e. assuming that all sink particles have zero age."
        )
        print(msg)
        sinks.birth_time = 0 | units.Myr


    print("Assigning groups to the sink particles...")
    sinks.initialised = False
    sinks = sinks.sorted_by_attribute('mass').reversed()
    Nsinks = len(sinks)
    for i, sink in enumerate(sinks):
        if i % int(Nsinks/5) == 0:
            print(f"Progress: {100*i/Nsinks:.2f}%")

        sink = assign_sink_group(
            sink,
            sinks,
            group_radius=grouping_distance,
            group_speed=grouping_speed,
            group_age=grouping_age,
        )

    print(f"{len(set(sinks.in_group))} groups found.")

    # Sanity check: each sink particle must be in a group.
    ungrouped_sinks = sinks.select_array(
        lambda x: x <= 0, ['in_group']
    )
    if not ungrouped_sinks.is_empty():
        print(
            "ERROR: There exist ungrouped sinks! Check group assignment."
        )
        exit()

    print("Forming stars from sink groups...")
    Ngroup = sinks.in_group.max()
    stars = Particles()
    for i in range(Ngroup):
        if i % int(Ngroup/5) == 0:
            print(f"Progress: {100*i/Ngroup:.2f}%")

        i += 1
        new_stars = form_stars_from_group(
            group_index=i,
            sink_particles=sinks,
            lower_mass_limit=mass_lower,
            upper_mass_limit=mass_upper,
            randomseed=numpy.random.randint(2**32-1),
            shrink_sinks=False
        )
        if new_stars is not None:
            stars.add_particles(new_stars)

    Nstars = len(stars)
    Mstars = stars.total_mass().in_(units.MSun)
    print(f"Number of stars: {Nstars}")
    print(f"Total star mass: {Mstars}")
    print(f"Most massive star: {stars.mass.max().in_(units.MSun)}")
    print(f"Average star mass: {Mstars/Nstars}")
    write_set_to_file(stars, starfilename, 'amuse')

    print(f"Star file written as {starfilename}")


if __name__ == '__main__':
    main()
