#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

"""
Module genAlg.py provides a genetic algorithm for optimization.

Copyright (c) 2018 Alan D. Freed

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

# Module metadata
__version__ = "0.1.1"
__date__ = "2017-09-20"
__author__ = "Alan D. Freed"
__author_email__ = "adfreed@tamu.edu"


"""
    Specie is a parent class to be inherited by all applications where a model
    is to be parameterized by a genetic algorithm.  Here is where the user
    creates their 'model' for whose 'parameters' a user seeks to quantify
    through this optimization scheme.

    The model is to assign its parameters array so it indexes from 1.  This is
    because index 0 is reserved for internal use by the optimizer to determine
    its optimal objective function.  This is a vital feature of this scheme.

    Constructor

        s = Specie(expTypes, ptsInExp, ctrlsInExp, respsInExp)

        where

            expTypes        number of experiment types with separate data sets
            ptsInExp        number of datum points in experiment[index]
            ctrlsInExp      number of controlled variables in experiment[index]
            respsInExp      number of  response  variables in experiment[index]

        with ptsInExp, ctrlsInExp, respsInExp being integer arrays.

    Attributes

    These attributes are accessed by the genetic algorithm.

        nExp    integer:       number of different experiment types
        nPts    integer array: number of data points in each experiment type
        nCtl    integer array: number of control variables per experiment type
        nRsp    integer array: number of response variables per experiment type
        expC    array of matrix: sized at exp[i] X [ctrl[i] X atPt[i]]
                                 i indexes from 0 to nExp-1
                                 array is of type object (an array of matrices)
                                 matrix [nCtl[i] X nPts[i]] is of type float
        expR    array of matrix: sized at exp[i] X [resp[i] X atPt[i]]
                                 i indexes from 0 to nExp-1
                                 array is of type object (an array of matrices)
                                 matrix [nRsp[i] X nPts[i]] is of type float
        modR    array of matrix: sized at exp[i] X [resp[i] X atPt[i]]
                                 i indexes from 0 to nExp-1
                                 array is of type object (an array of matrices)
                                 matrix [nRsp[i] X nPts[i]] is of type float

    Methods

    The following methods are called prior to calling the genetic algorithm.

        data = s.newExpCtrlData()
            data            creates a new array of matrices to be filled in by
                            a user, e.g., data[i][j,k], where:
                            i)  array indexes over the number of experiments
                            ii) matrix has rows that index over the number
                                of controlled variables for the specified
                                experiment, and columns that index over the
                                number of data points for that experiment
        s.assignExpCtrlData(data)
            data            assigns a filled in data structure to self.expC
        data = s.newExpRespData()
            data            creates a new array of matrices to be filled in by
                            a user, e.g., data[i][j,k], where:
                            i)  array indexes over the number of experiments
                            ii) matrix has rows that index over the number
                                of response variables for the specified
                                experiment, and columns that index over the
                                number of data points for that experiment
        s.assignExpRespData(data)
            data            assigns a filled in data structure to self.expR

    The following method is called from within the genetic algorithm.

        s.model(parameters) solves the user's model using the parameters given
                            and places the response results at self.modR
            parameters      array of model parameters managed by the genetic
                            algorithm: start indexing your model's parameters
                            with index 1, as index 0 is reserved for the
                            parameter within the objective function itself
"""


class Specie(object):

    def __init__(self, expTypes, ptsInExp, ctrlsInExp, respsInExp):
        # Must be called via super by all child objects.

        # The number of experiments or experiment types to be analyzed.
        # Each experiment/type may be described with a different model.
        if expTypes > 0:
            self.nExp = int(expTypes)
        else:
            raise RuntimeError("Number of 'expTypes' must be positive.")

        # Supplies number of datum points that comprise each experiment type.
        if isinstance(ptsInExp, np.ndarray):
            if ptsInExp.size == expTypes:
                self.nPts = np.zeros(self.nExp, dtype=int)
            else:
                raise RuntimeError(
                          "The ength of 'ptsPerExp' did not equal 'expTypes'.")
        else:
            raise RuntimeError("Argument 'ptsPerExp' must be a numpy array.")
        for i in range(self.nExp):
            self.nPts[i] = int(ptsInExp[i])

        # Supplies number of controlled variables for each experiment type.
        if isinstance(ctrlsInExp, np.ndarray):
            if ctrlsInExp.size == expTypes:
                self.nCtl = np.zeros(self.nExp, dtype=int)
            else:
                raise RuntimeError(
                     "The ength of 'ctrlsPerExp' did not equal 'experiments'.")
        else:
            raise RuntimeError("Argument 'ctrlsPerExp' must be a numpy array.")
        for i in range(self.nExp):
            self.nCtl[i] = int(ctrlsInExp[i])

        # Supplies the number of response variables for each experiment.
        if isinstance(respsInExp, np.ndarray):
            if respsInExp.size == self.nExp:
                self.nRsp = np.zeros(self.nExp, dtype=int)
            else:
                raise RuntimeError(
                     "The ength of 'respsPerExp' did not equal 'experiments'.")
        else:
            raise RuntimeError("Argument 'respsPerExp' must be a numpy array.")
        for i in range(self.nExp):
            self.nRsp[i] = int(respsInExp[i])

        return  # new object: this is the constructor

    # Create data structure for filling in the experimental control variables.
    def newExpCtrlData(self):
        # Create the data structures.
        control = np.empty(self.nExp, dtype=object)
        for i in range(self.nExp):
            control[i] = np.zeros(shape=(self.nCtl[i], self.nPts[i]),
                                  dtype=float)
        return control

    # Assign filled-in data structure for the experimental control variables.
    def assignExpCtrlData(self, controlData):
        self.expC = controlData

    # Create data structure for filling in the experimental response variables.
    def newExpRespData(self):
        # Create the data structures.
        response = np.empty(self.nExp, dtype=object)
        for i in range(self.nExp):
            response[i] = np.zeros(shape=(self.nRsp[i], self.nPts[i]),
                                   dtype=float)
        return response

    # Assign filled-in data structure for the experimental response variables.
    def assignExpRespData(self, responseData):
        self.expR = responseData

    # The model's interface to the genetic algorithm.
    def model(self, parameters):
        # Must be called via super by all child objecets.

        if not isinstance(parameters, np.ndarray):
            raise RuntimeError("Sent 'parameters' must be a numpy array.")

        #  Create the response data structure.
        self.modR = np.empty(self.nExp, dtype=object)
        for i in range(self.nExp):
            self.modR[i] = np.zeros(shape=(self.nRsp[i], self.nPts[i]),
                                    dtype=float)
        #  Return control back to the child object to fill in these data.
        return  # nothing here
