# Jack C. Cook
# Friday, December 10, 2021

import copy
import numpy as np
import cv2
import ghedt as dt


class Constraints:
    def __init__(self, length: float = None, width: float = None,
                 B: float = None, B_min: float = None,
                 B_max_x: float = None, B_max_y: float = None,
                 outer_constraints: list = None, no_go: list = None):
        # Spacing parameters in meters
        self.B = B
        self.B_max_x = B_max_x
        self.B_max_y = B_max_y
        self.B_min = B_min
        # Length and width constraints
        self.length = length
        self.width = width
        # Outer constraints described as a polygon
        self.outer_constraints = outer_constraints
        # TODO: Handling for a list or a list of lists to occur later
        # Note: the entirety of the no-go zone should fall inside the
        # outer_constraints
        self.no_go = no_go

    def check_inputs(self, method):
        # The required instances for the near-square design is self.B
        if method == 'near-square':
            assert self.B is not None
            assert self.length is not None
        elif method == 'rectangle':
            assert self.width is not None
            assert self.length is not None
            assert self.B_min is not None
            assert self.B_max_x is not None
        elif method == 'bi-rectangle' or 'bi-zoned':
            assert self.width is not None
            assert self.length is not None
            assert self.B_min is not None
            assert self.B_max_x is not None
            assert self.B_max_y is not None

        return


def remove_cutout(coordinates, boundary=None, remove_inside=True,
                  keep_contour=True):
    if boundary is None:
        boundary = []

    # cv2.pointPolygonTest only takes integers, so we scale by 10000 and then
    # scale back to keep precision
    scale = 10000.
    coordinates = dt.coordinates.scale_coordinates(coordinates, scale)
    boundary = dt.coordinates.scale_coordinates(boundary, scale)

    _boundary = np.array(boundary, dtype=np.uint64)

    # https://stackoverflow.com/a/50670359/11637415
    # Positive - point is inside the contour
    # Negative - point is outside the contour
    # Zero - point is on the contour

    inside_points_idx = []
    outside_points_idx = []
    boundary_points_idx = []
    for i in range(len(coordinates)):
        coordinate = coordinates[i]
        dist = cv2.pointPolygonTest(_boundary, coordinate, False)
        if dist > 0.0:
            inside_points_idx.append(i)
        elif dist < 0.0:
            outside_points_idx.append(i)
        elif dist == 0.0:
            boundary_points_idx.append(i)

    new_coordinates = []
    for i in range(len(coordinates)):
        # if we want to remove inside points and keep contour points
        if remove_inside and keep_contour:
            if i in inside_points_idx:
                continue
            else:
                new_coordinates.append(coordinates[i])
        # if we want to remove inside points and remove contour points
        elif remove_inside and not keep_contour:
            if i in inside_points_idx or i in boundary_points_idx:
                continue
            else:
                new_coordinates.append(coordinates[i])
        # if we want to keep outside points and remove contour points
        elif not remove_inside and not keep_contour:
            if i in outside_points_idx or i in boundary_points_idx:
                continue
            else:
                new_coordinates.append(coordinates[i])
        # if we want to keep outside points and keep contour points
        else:
            if i in outside_points_idx:
                continue
            else:
                new_coordinates.append(coordinates[i])

    new_coordinates = dt.coordinates.scale_coordinates(new_coordinates, 1/scale)

    return new_coordinates


def determine_largest_rectangle(property_boundary):
    x_max = 0
    y_max = 0
    for i in range(len(property_boundary)):
        x, y = property_boundary[i]
        if x > x_max:
            x_max = copy.deepcopy(x)
        if y > y_max:
            y_max = copy.deepcopy(y)

    rectangle = [[0, 0],
                 [x_max, 0],
                 [x_max, y_max],
                 [0, y_max],
                 [0, 0]]

    return rectangle
