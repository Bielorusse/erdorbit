"""
Graphics module for erdorbit.
"""

# third party imports
import numpy as np
from shapely.geometry import LineString


def resize_drawing_to_fit_canvas(input_coords, canvas_height, drawing_size_factor):
    """
    This function resizes the drawing to fit the canvas size
    Input:
        input_coords            np array of floats shape (:, 3)
            x, y, z, position in km
        canvas_height           float
        drawing_size_factor     float
    Output:
        np array of floats shape (:, 3)
            x, y, z, position in canvas dimension
    """

    # assuming size of object to draw is 2 times largest value in positions array
    drawing_size = np.max(np.abs(input_coords)) * 2

    return input_coords * canvas_height / drawing_size * drawing_size_factor


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


def get_projected_coords(
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
    Converts the 3D positions of the orbit to projected 2D positions intended for drawing.
    Input:
        pos             np array of floats shape (:, 3)
            array of 3D cartesian positions of the orbiting object
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
    Output:
        positions_2d    np array of floats shape (:, 2)
            array of 2D cartesian positions of the orbiting object
    """

    positions_2d = planar_projection(positions, alpha, beta)

    positions_2d = rotate_coordinates(positions_2d, delta)

    # translate coordinates
    positions_2d[:, 0] += x_translation
    positions_2d[:, 1] += y_translation

    return positions_2d


def from_cartesian_to_computer_coords(positions_2d, canvas_width, canvas_height):
    """
    Conversion from cartesian to computer coordinates.
    Input:
        positions_2d    np array of floats shape (:, 2)
            array of 2D cartesian positions of the orbiting object
        canvas_width    float
        canvas_height   float
    Output:
        positions_2d    np array of floats shape (:, 2)
            array of 2D computer positions of the orbiting object
    """

    # conversion from cartesian to computer coordinates.
    positions_2d[:, 0] += canvas_width / 2
    positions_2d[:, 1] = -positions_2d[:, 1] + canvas_height / 2

    return positions_2d


def write_csv_file(positions_2d, filename):
    """
    Write list of 2d positions in csv file.
    Input:
        -positions_2d   np array of floats shape (:, 2)
        -filename       str
    """
    with open(filename, "w") as outfile:
        for i in range(positions_2d.shape[0]):
            outfile.write("{}, {}\n".format(positions_2d[i, 0], positions_2d[i, 1]))


def write_svg_file(
    positions_2d,
    filename,
    canvas_width,
    canvas_height,
    title=None,
    description=None,
):
    """
    Write list of 2d positions in svg file as lines.
    Input:
        -positions_2d   np array of floats shape (:, 2)
        -filename       str
        -canvas_width   int
        -canvas_height  int
        -title          str
        -description    str
    """
    with open(filename, "w") as outfile:

        # write headers
        outfile.write("<?xml version='1.0' encoding='utf-8'?>\n")
        outfile.write(
            "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' width='{}' height='{}'>\n".format(
                canvas_width, canvas_height
            )
        )
        if title:
            outfile.write("<title>{}</title>\n".format(title))
        if description:
            outfile.write("<description>{}</description>\n".format(description))

        # loop through 2d positions
        for i in range(positions_2d.shape[0] - 1):

            # create linestring
            line = LineString([positions_2d[i, :], positions_2d[i + 1, :]])
            outfile.write("{}\n".format(line.svg()))

        # write footer
        outfile.write("</svg>")
