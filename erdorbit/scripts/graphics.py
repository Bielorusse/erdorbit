"""
Graphics module for erdorbit.
"""

# third party imports
import numpy as np
from processing_py import *
from tqdm import tqdm


def resize_drawing_to_fit_canvas(positions, canvas_height, DRAWING_SIZE_FACTOR):
    """
    This function resizes the drawing to fit the canvas size
    Input:
        input_coords            np array of floats shape (:, 3)
            x, y, z, position in km
        canvas_height           float
        DRAWING_SIZE_FACTOR     float
    Output:
        output_coords           np array of floats shape (:, 3)
            x, y, z, position in canvas dimension
    """

    # assuming size of object to draw is 2 times largest value in positions array
    drawing_size = np.max(np.abs(positions)) * 2

    return positions * canvas_height / drawing_size * DRAWING_SIZE_FACTOR


def planar_projection(input_array, alpha, beta):
    """
    Projecting an array of 3D coordinates on a plane
    Inputs:
        input_array     np array of floats shape (:, 3)
            3D position expressed as x, y, z
        alpha           float
            rotation of the projection plane around intersection of itself and Oxy plane (radians)
        beta            float
            rotation of the projection plane around Oz (radians)
    Outputs:
        output_array    np array of floats shape (:, 2)
            2D projection of input positions expressed as x, y
    """

    output_array = np.zeros((input_array.shape[0], 2))
    for i in range(input_array.shape[0]):
        output_array[i, 0] = (
            np.cos(beta) * input_array[i, 0] - np.sin(beta) * input_array[i, 1]
        )
        output_array[i, 1] = (
            -np.sin(beta) * np.sin(alpha) * input_array[i, 0]
            - np.cos(beta) * np.sin(alpha) * input_array[i, 1]
            + np.cos(alpha) * input_array[i, 2]
        )

    return output_array


def rotate_coordinates(input_coords, delta):
    """
    This functions rotates 2D coordinates around the center of the reference
    frame from a 'delta' angle.
    Input:
        input_coords    np array of floats shape (:, 2)
            2D coordinates expressed as x, y
        delta           float
            rotation of the projected vector around the normal to the projection plane, in radians
    Output:
        output_coords   np array of floats shape (:, 2)
            2D coordinates expressed as x, y
    """

    output_coords = np.zeros((input_coords.shape[0], 2))
    for i in range(input_coords.shape[0]):
        output_coords[i, 0] = (
            np.cos(delta) * input_coords[i][0] + np.sin(delta) * input_coords[i][1]
        )
        output_coords[i, 1] = (
            np.cos(delta) * input_coords[i][1] - np.sin(delta) * input_coords[i][0]
        )

    return output_coords


def draw_orbit(
    positions,
    alpha,
    beta,
    delta,
    x_translation,
    y_translation,
    canvas_width,
    canvas_height,
):
    """
    This function draws the trajectory of an orbiting object on the canvas
    This function calls several previously defined functions
    Input:
        pos             np array of floats shape (:, 3)
            array of 3D cartesian positions the orbiting object
        alpha           float
            projection plane angle
        beta            float
            projection plane angle
        delta           float
            drawing rotation angle
        x_translation   float
            drawing translation along x-axis
        y_translation   float
            drawing translation along y-axis
        canvas_width    float
        canvas_height   float
    """

    app = App(int(canvas_width), int(canvas_height))  # create window: width, height
    app.background(255, 255, 255)  # set white background
    app.stroke(0, 0, 0)  # set stroke color to black

    positions_2d = planar_projection(positions, alpha, beta)

    positions_2d = rotate_coordinates(positions_2d, delta)

    # translate coordinates
    positions_2d[:, 0] += x_translation
    positions_2d[:, 1] += y_translation

    # conversion from cartesian to computer coordinates.
    positions_2d[:, 0] += canvas_width / 2
    positions_2d[:, 1] = -positions_2d[:, 1] + canvas_height / 2

    # drawing line
    for i in tqdm(range(len(positions_2d) - 1)):
        app.line(
            positions_2d[i][0],
            positions_2d[i][1],
            positions_2d[i + 1][0],
            positions_2d[i + 1][1],
        )

    app.redraw()  # refresh drawing
