"""
The purpose of this script is to draw erdorbit with the AxiDraw plotter.
"""

# standard imports
import os
import sys
import argparse
import configparser

# third party imports
import numpy as np

# local imports
from computation import compute_positions
from graphics import resize_drawing_to_fit_canvas
from graphics import get_projected_coords
from graphics import from_cartesian_to_computer_coords
from graphics import write_csv_file
from graphics import write_svg_file
from drawing import draw_orbit


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
    ROT_PLANET = (
        float(config["calculation"]["ROT_PLANET"]) / 86164.1004
    )  # translation from deg/sidereal-day to deg/s
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
    alpha = float(config["drawing"]["alpha"]) * np.pi / 180
    beta = float(config["drawing"]["beta"]) * np.pi / 180
    delta = float(config["drawing"]["delta"]) * np.pi / 180
    x_translation = float(config["drawing"]["x_translation"])
    y_translation = float(config["drawing"]["y_translation"])
    canvas_height = float(config["drawing"]["canvas_height"])
    canvas_width = float(config["drawing"]["canvas_width"])
    DRAWING_SIZE_FACTOR = float(config["drawing"]["DRAWING_SIZE_FACTOR"])
    csv_output_file = config["drawing"]["csv_output_file"]
    svg_output_file = config["drawing"]["svg_output_file"]
    display_p5_preview = config["drawing"]["display_p5_preview"]

    # computing positions
    positions = compute_positions(
        a, e, inc, RAAN, om, MU_PLANET, ROT_PLANET, number_of_steps, step_time
    )

    # resizing to canvas
    positions = resize_drawing_to_fit_canvas(
        positions, canvas_height, DRAWING_SIZE_FACTOR
    )

    # get 2d projected coords
    positions_2d = get_projected_coords(
        positions,
        alpha,
        beta,
        delta,
        x_translation,
        y_translation,
        canvas_width,
        canvas_height,
    )

    # coordinates conversion from cartesian to computer
    positions_2d = from_cartesian_to_computer_coords(
        positions_2d, canvas_width, canvas_height
    )

    if csv_output_file != "None":
        write_csv_file(positions_2d, csv_output_file)

    if svg_output_file != "None":
        write_svg_file(positions_2d, svg_output_file, canvas_width, canvas_height)

    if display_p5_preview.lower() in ["true", "t", "y", "yes"]:
        draw_orbit(positions_2d, canvas_width, canvas_height)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config")
    args = parser.parse_args()

    erdorbit(args.config)
