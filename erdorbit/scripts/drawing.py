"""
Module for drawing 2d projected coordinates intended for erdorbit.
"""

# third party imports
from processing_py import *
from tqdm import tqdm


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
