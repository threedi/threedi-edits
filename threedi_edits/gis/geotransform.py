# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 14:51:40 2021

@author: chris.kerklaan
"""
import math


class GeoTransform(tuple):
    def shifted(self, geometry, inflate=False):
        """
        Return shifted geo transform.

        :param geometry: geometry to match
        :param inflate: inflate to nearest top-left grid intersection.
        """
        values = list(self)
        index = self.get_indices(geometry, inflate=inflate)[1::-1]
        values[0], values[3] = self.get_coordinates(index)
        return self.__class__(values)

    def scaled(self, w, h):
        """
        Return shifted geo transform.

        :param f: scale the cellsize by this factor
        """
        p, a, b, q, c, d = self
        return self.__class__([p, a * w, b * h, q, c * w, d * h])

    def get_coordinates(self, indices):
        """Return x, y coordinates.

        :param indices: i, j tuple of integers or arrays.

        i corresponds to the y direction in a non-skew grid.
        """
        p, a, b, q, c, d = self
        i, j = indices
        return p + a * j + b * i, q + c * j + d * i

    def get_indices(self, geometry, inflate=False):
        """
        Return array indices tuple for geometry.

        :param geometry: geometry to subselect
        :param inflate: inflate envelope to grid, to make sure that
            the entire geometry is contained in resulting indices.
        """
        # spatial coordinates
        x1, x2, y1, y2 = geometry.GetEnvelope()

        # inverse transformation
        p, a, b, q, c, d = self
        e, f, g, h = get_inverse(a, b, c, d)

        # apply to envelope corners
        f_lo, f_hi = (math.floor, math.ceil) if inflate else (round, round)

        X1 = int(f_lo(e * (x1 - p) + f * (y2 - q)))
        Y1 = int(f_lo(g * (x1 - p) + h * (y2 - q)))
        X2 = int(f_hi(e * (x2 - p) + f * (y1 - q)))
        Y2 = int(f_hi(g * (x2 - p) + h * (y1 - q)))

        # prevent zero dimensions in case of inflate
        if inflate:
            if X1 == X2:
                X2 += 1
            if Y1 == Y2:
                Y1 -= 1

        return X1, Y1, X2, Y2

    def get_slices(self, geometry):
        """
        Return array slices tuple for geometry.

        :param geometry: geometry to subselect
        """
        x1, y1, x2, y2 = self.get_indices(geometry)
        return slice(y1, y2), slice(x1, x2)

    def get_window(self, geometry):
        """
        Return window dictionary for a geometry.

        :param geometry: geometry to subselect
        """
        x1, y1, x2, y2 = self.get_indices(geometry)
        return {"xoff": x1, "yoff": y1, "xsize": x2 - x1, "ysize": y2 - y1}


def get_inverse(a, b, c, d):
    """Return inverse for a 2 x 2 matrix with elements (a, b), (c, d)."""
    D = 1 / (a * d - b * c)
    return d * D, -b * D, -c * D, a * D
