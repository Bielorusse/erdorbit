"""
The purpose of this script is to draw erdorbit with the AxiDraw plotter.
"""

# standard imports
import os
import sys
import argparse
import configparser

# local imports
from computation import compute_positions
from graphics import resize_drawing_to_fit_canvas
from graphics import draw_orbit


def erdorbit(config_file=None):
    """
    Draw erdorbit with the AxiDraw plotter.
    Input:
        -config_file    str
    """

    # read config
    config = configparser.ConfigParser()
    default_config_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "default_config.ini"
    )
    config.read(default_config_file)  # first read default config
    if not config_file is None:  # replace with new parameters if custom config given
        config.read(config_file)

    # calculation parameters
    MU_PLANET = float(config["calculation"]["MU_PLANET"])
    ROT_PLANET = float(config["calculation"]["ROT_PLANET"])
    duration = (
        float(config["calculation"]["duration"]) * 86400
    )  # conversion from days to seconds
    if "step_time" in config["calculation"]:
        step_time = float(config["calculation"]["step_time"])
        number_of_steps = duration / step_time
    elif "number_of_steps" in config["calculation"]:
        number_of_steps = int(config["calculation"]["number_of_steps"])
        step_time = duration / number_of_steps
    else:
        sys.exit("Error reading config")
    a = float(config["calculation"]["a"])
    e = float(config["calculation"]["e"])
    inc = float(config["calculation"]["i"])
    RAAN = float(config["calculation"]["RAAN"])
    om = float(config["calculation"]["om"])

    # drawing parameters
    alpha = float(config["drawing"]["alpha"])
    beta = float(config["drawing"]["beta"])
    delta = float(config["drawing"]["delta"])
    x_translation = float(config["drawing"]["x_translation"])
    y_translation = float(config["drawing"]["y_translation"])
    canvas_height = float(config["drawing"]["canvas_height"])
    canvas_width = float(config["drawing"]["canvas_width"])
    DRAWING_SIZE_FACTOR = float(config["drawing"]["DRAWING_SIZE_FACTOR"])

    # computing positions
    positions = compute_positions(
        a, e, inc, RAAN, om, MU_PLANET, ROT_PLANET, number_of_steps, step_time
    )

    # resizing to canvas
    positions = resize_drawing_to_fit_canvas(
        positions, canvas_height, DRAWING_SIZE_FACTOR
    )

    # drawing orbit
    draw_orbit(
        positions,
        alpha,
        beta,
        delta,
        x_translation,
        y_translation,
        canvas_width,
        canvas_height,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config")
    args = parser.parse_args()

    erdorbit(args.config)
