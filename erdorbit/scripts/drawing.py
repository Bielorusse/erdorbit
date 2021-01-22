"""
Module for drawing 2d projected coordinates intended for erdorbit.
"""

# standard imports
import argparse

# third party imports
from processing_py import *
from tqdm import tqdm
import numpy as np


def draw_orbit(positions_2d, canvas_width, canvas_height):
    """
    Draw projected 2d coordinates on canvas.
    Input:
        positions_2d    np array of floats shape (:, 2)
            array of 2D computer positions the orbiting object
        canvas_width    int
        canvas_height   int
    """

    # creating some constants for drawing with processing
    app = App(int(canvas_width), int(canvas_height))  # create window: width, height
    app.background(255, 255, 255)  # set white background
    app.stroke(0, 0, 0)  # set stroke color to black

    # drawing line
    for i in tqdm(range(len(positions_2d) - 1)):
        app.line(
            positions_2d[i][0],
            positions_2d[i][1],
            positions_2d[i + 1][0],
            positions_2d[i + 1][1],
        )

    app.redraw()  # refresh drawing


def read_csv(input_file):
    """
    Read 2d positions csv file and store contents in a numpy array.
    Input:
        -input_file     str
    """
    positions = []  # initiate positions list
    with open(input_file, "r") as infile:
        for row in infile:
            values = [float(col.strip()) for col in row.split(",")]  # read values
            positions.append(values)  # store values in array
    return np.asarray(positions)  # return positions list as numpy array


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("-csv", "--input_csv_file", required=True)
    parser.add_argument("-wi", "--canvas_width")
    parser.add_argument("-he", "--canvas_height")
    args = parser.parse_args()

    # read csv
    positions_2d = read_csv(args.input_csv_file)

    # get canvas width and height
    if args.canvas_width is not None and args.canvas_height is not None:
        canvas_width = int(args.canvas_width)
        canvas_height = int(args.canvas_height)
    else:  # compute width and height based on 2d positions
        canvas_width = int(np.max(np.abs(positions_2d[:, 0])))
        canvas_height = int(np.max(np.abs(positions_2d[:, 1])))

    # draw orbit
    draw_orbit(positions_2d, canvas_width, canvas_height)
