"""
Computation module for erdorbit.
"""

# third party imports
import numpy as np


def compute_positions(
    a, e, inc, RAAN, om, MU_PLANET, ROT_PLANET, number_of_steps, step_time
):
    """
    This function computes the positions of the orbiting object in a ECEF reference frame
    Inputs:
        orbital parameters          [float, float, float, float, float]
            a           semi-major axis (km)
            e           eccentricity (unitless)
            inc         inclination (deg)
            RAAN        right ascension of the ascending node (deg)
            om          argument of periapsis (deg)
        planetary constants         [float, float]
            MU_PLANET   planet gravitational parameter (km3/s2)
            ROT_PLANET  planet rotational velocity (deg/s)
        number_of_steps             int
        step_time (s)               int
    Ouputs:
        positions (x, y, z, in km)  [[float, float, float], ...]
    """

    positions = []

    # computing positions for each simulation step
    for i in range(number_of_steps):

        time_elapsed = i * step_time

        # computing this position vector for this step
        positions.append(
            from_orbital_to_cartesian_coordinates(
                a, e, inc, RAAN, om, time_elapsed, MU_PLANET
            )
        )

        # rotating frame around z axis for conversion from J2000 to ECEF reference frame
        positions[-1] = rotate_frame_around_z(positions[-1], ROT_PLANET * time_elapsed)

    return np.asarray(positions)


def rotate_frame_around_z(input_vector, angle):
    """
    Converts coordinates to rotated reference frame.
    Inputs:
        input_vector    [float, float, float]
        angle (deg)     float
    Outputs:
        output_vector   [float, float, float]
    """

    angle = angle * np.pi / 180

    output_vector = [
        np.cos(angle) * input_vector[0] + np.sin(angle) * input_vector[1],
        np.cos(angle) * input_vector[1] - np.sin(angle) * input_vector[0],
        input_vector[2],
    ]

    return output_vector


def from_orbital_to_cartesian_coordinates(a, e, inc, raan, om, t, mu):
    """
    Converting from orbital parameters to cartesian coordinates.
    Input:
        a       float                   semi-major axis (km)
        e       float                   eccentricity (-)
        inc     float                   inclination (deg)
        raan    float                   right ascension of the ascending node (deg)
        om      float                   argument of periapsis (deg)
        t       float                   time spent since passage at periapsis (s)
        mu      float                   gravitational parameter of the central body (km3/s2)
    Outputs:
        pos     [float, float, float]   x, y, z, in km
    """

    # converting angles from degrees to radians
    inc = inc * np.pi / 180
    raan = raan * np.pi / 180
    om = om * np.pi / 180

    # computing mean anomaly
    n = np.sqrt(mu / np.power(a, 3.0))
    M = n * t

    # computing eccentric anomaly
    E = [M]
    for j in range(100):
        E.append(E[j] + (M - E[j] + e * np.sin(E[j])) / (1 - e * np.cos(E[j])))
        if abs(E[j + 1] - E[j]) < 1e-8:
            E = E[j + 1]
            break

    # computing true anomaly
    nu = (
        2
        * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
        % (np.pi * 2)
    )

    # computing radius
    r = a * (1 - np.power(e, 2.0)) / (1 + e * np.cos(nu))

    # computing position vector
    pos = np.asarray(
        (
            r
            * (
                np.cos(om + nu) * np.sin(raan)
                + np.sin(om + nu) * np.cos(raan) * np.cos(inc)
            ),
            r * (np.sin(om + nu) * np.sin(inc)),
            r
            * (
                np.cos(om + nu) * np.cos(raan)
                - np.sin(om + nu) * np.sin(raan) * np.cos(inc)
            ),
        )
    )

    return pos
