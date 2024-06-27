#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import numpy as np
from genAlgs import GenAlg
from species import Specie
# for creating graphics
from matplotlib import pyplot as plt
from matplotlib import rc

"""
A test code written for genAlg.py.
"""

#   The model parameters for the Garofalo/Miller steady-state creep model are:
#       P[0] = p      #  the exponent to computed error for objective function
#       P[1] = Tm     #  melting temperature in K, which does not vary
#       P[2] = Q      #  activation energy for self diffusion
#       P[3] = Tt/Tm  #  location of the transition temperature
#       P[4] = A      #  power-law coefficient
#       P[5] = C      #  power-law breakdown stress
#       P[6] = n      #  power-law exponent


#   Garafolo/Miller - Zener/Holloman theory for uniaxial steady-state creep.
class ssCreep(Specie):

    def __init__(self, expTypes, ptsInExp, ctrlsInExp, respsInExp):
        # create an instance of the parent class
        super(ssCreep, self).__init__(expTypes, ptsInExp,
                                      ctrlsInExp, respsInExp)
        # create other attributes
        self.sqrt3 = math.sqrt(3.0)
        return  # nothing

    # Create data structure for filling in the experimental control variables.
    def newExpCtrlData(self):
        # Create the data structures.
        control = super(ssCreep, self).newExpCtrlData()
        return control

    # Place data structure for the experimental control variables.
    def assignExpCtrlData(self, controlData):
        super(ssCreep, self).assignExpCtrlData(controlData)

    # Create data structure for filling in the experimental response variables.
    def newExpRespData(self):
        # Create the data structures.
        response = super(ssCreep, self).newExpRespData()
        return response

    # Place data structure for the experimental response variables.
    def assignExpRespData(self, responseData):
        super(ssCreep, self).assignExpRespData(responseData)

    def theta(self, P, T):
        #  universal gas constant
        R = 8.314
        # physical parameter
        Tm = P[1]
        # diffusion related parameters
        Q = P[2]
        Tt = P[3] * Tm
        # Miller's function for thermal activation
        if T > Tt:
            fn = math.exp(-Q / (R * T))
        else:
            fn = math.exp((-Q / (R * Tt)) * (math.log(Tt / T) + 1.0))
        return fn

    def Zener(self, P, S):
        # creep model parameters
        A = 10**(P[4])  # these data are logged during optimization
        C = P[5]
        n = P[6]
        # Garafolo's creep law
        fn = A * (math.sinh(S / (self.sqrt3 * C))) ** n
        return fn

    # Zener-Hollaman model for steady-state creep in uniaxial extension
    def model(self, parameters):
        # create the data structue that holds the response variables
        super(ssCreep, self).model(parameters)
        # fill in this data structure

        for i in range(self.nExp):
            self.modR[i] = np.zeros(shape=(self.nRsp[i], self.nPts[i]),
                                    dtype=float)

        for i in range(self.nExp):
            for j in range(self.nRsp[i]):
                for k in range(self.nPts[i]):
                    # the controlled variables are
                    T = self.expC[i][0, k]   # temperature in K
                    S = self.expC[i][1, k]   # stress in MPa
                    # the model response is
                    creepRate = (self.theta(parameters, T) *
                                 self.Zener(parameters, S) / self.sqrt3)
                    self.modR[i][j, k] = creepRate
        return  # nothing


def runGenAlg():
    # temperatures in K at which creep tests were run
    temp = ([1089.0, 1255.0, 1422.0, 1589.0])

    # input the creep data as tuples:
    # (stress (MPa), steady-state creep rate (1/s))
    creepData1 = []  # at 1089 K
    creepData2 = []  # at 1255 K
    creepData3 = []  # at 1422 K
    creepData4 = []  # at 1589 K
    # data at 1089 K: 13 data pairs
    creepData1.append((50.1, 0.37e-6))
    creepData1.append((61.2, 0.94e-6))
    creepData1.append((66.9, 1.7e-6))
    creepData1.append((72.5, 3.5e-6))
    creepData1.append((78.1, 5.8e-6))
    creepData1.append((84.6, 8.0e-6))
    creepData1.append((92.5, 15.0e-6))
    creepData1.append((100.5, 31.0e-6))
    creepData1.append((54.9, 0.79e-6))
    creepData1.append((58.7, 1.5e-6))
    creepData1.append((63.1, 2.5e-6))
    creepData1.append((72.1, 6.6e-6))
    creepData1.append((82.7, 22.0e-6))
    # data at 1255 K: 14 data pairs
    creepData2.append((20.0, 0.34e-6))
    creepData2.append((21.2, 1.3e-6))
    creepData2.append((24.1, 1.7e-6))
    creepData2.append((27.4, 2.9e-6))
    creepData2.append((30.8, 4.8e-6))
    creepData2.append((34.5, 9.8e-6))
    creepData2.append((38.3, 22.0e-6))
    creepData2.append((42.6, 51.0e-6))
    creepData2.append((25.6, 0.73e-6))
    creepData2.append((30.1, 4.0e-6))
    creepData2.append((33.1, 6.2e-6))
    creepData2.append((36.5, 11.0e-6))
    creepData2.append((40.7, 24.0e-6))
    creepData2.append((45.4, 82.0e-6))
    # data at 1422 K: 11 data pairs
    creepData3.append((7.93, 0.3e-6))
    creepData3.append((11.5, 1.7e-6))
    creepData3.append((14.5, 4.8e-6))
    creepData3.append((18.6, 14.0e-6))
    creepData3.append((10.5, 2.0e-6))
    creepData3.append((11.1, 2.2e-6))
    creepData3.append((11.8, 2.4e-6))
    creepData3.append((12.9, 2.8e-6))
    creepData3.append((13.9, 3.8e-6))
    creepData3.append((15.4, 9.1e-6))
    creepData3.append((18.1, 21.0e-6))
    # data at 1589 K: 11 data pairs
    creepData4.append((3.7, 0.35e-6))
    creepData4.append((4.13, 0.52e-6))
    creepData4.append((5.01, 1.5e-6))
    creepData4.append((5.92, 3.3e-6))
    creepData4.append((7.45, 7.3e-6))
    creepData4.append((9.38, 20.0e-6))
    creepData4.append((11.6, 65.0e-6))
    creepData4.append((5.04, 2.8e-6))
    creepData4.append((5.67, 4.2e-6))
    creepData4.append((6.61, 7.8e-6))
    creepData4.append((7.72, 10.0e-6))

    # create data arrays by temperature for plotting

    stress1 = np.array([x[0] for x in creepData1])
    creepRate1 = np.array([x[1] for x in creepData1])
    stress2 = np.array([x[0] for x in creepData2])
    creepRate2 = np.array([x[1] for x in creepData2])
    stress3 = np.array([x[0] for x in creepData3])
    creepRate3 = np.array([x[1] for x in creepData3])
    stress4 = np.array([x[0] for x in creepData4])
    creepRate4 = np.array([x[1] for x in creepData4])

    # establish the model to be parameterized

    expTypes = 4  # one experiment 'type' per temperature
    ctrlsInExp = np.zeros(expTypes, dtype=int)
    respsInExp = np.zeros(expTypes, dtype=int)
    ptsInExp = np.zeros(expTypes, dtype=int)
    for i in range(expTypes):
        ctrlsInExp[i] = 2
        respsInExp[i] = 1
    ptsInExp[0] = stress1.size
    ptsInExp[1] = stress2.size
    ptsInExp[2] = stress3.size
    ptsInExp[3] = stress4.size

    ss = ssCreep(expTypes, ptsInExp, ctrlsInExp, respsInExp)

    # Create the experimental data structures.

    expCtrl = ss.newExpCtrlData()
    for k in range(ptsInExp[0]):
        expCtrl[0][0, k] = temp[0]
        expCtrl[0][1, k] = stress1[k]
    for k in range(ptsInExp[1]):
        expCtrl[1][0, k] = temp[1]
        expCtrl[1][1, k] = stress2[k]
    for k in range(ptsInExp[2]):
        expCtrl[2][0, k] = temp[2]
        expCtrl[2][1, k] = stress3[k]
    for k in range(ptsInExp[3]):
        expCtrl[3][0, k] = temp[3]
        expCtrl[3][1, k] = stress4[k]
    ss.assignExpCtrlData(expCtrl)

    expResp = ss.newExpRespData()
    for k in range(ptsInExp[0]):
        expResp[0][0, k] = creepRate1[k]
    for k in range(ptsInExp[1]):
        expResp[1][0, k] = creepRate2[k]
    for k in range(ptsInExp[2]):
        expResp[2][0, k] = creepRate3[k]
    for k in range(ptsInExp[3]):
        expResp[3][0, k] = creepRate4[k]
    ss.assignExpRespData(expResp)

    # run the genetic algorithm to determine optimal parameters

    sigFig = 6      # number of significant figures
    probM = 0.005   # probability of gene mutation at conception
    probX = 0.85    # probability of chromosome splitting at conception
    probI = 0.001   # probability of immigrant migrating into next generation

    nameP = np.array(['ObjFn p', 'Tm', 'Q', 'Tt/Tm', 'log10(A)', 'C', 'n'])
    varyP = np.array([True, False, True, True, True, True, True])
    fixdP = np.array([2173.0])
    alien = np.array([0.000105, 3.56e5, 0.489, 11.01, 35.08, 4.53])
    minP = np.array([0.00001, 2.0e5, 0.45, 10.5, 10.0, 3.0])
    maxP = np.array([0.001, 4.5e5, 0.85, 13.5, 75.0, 6.0])

    ga = GenAlg(ss, sigFig, probM, probX, probI,
                nameP, varyP, fixdP, alien, minP, maxP)
    text_file = open("genAlgOutput.txt", "w")
    text_file.write(ga.report())
    for i in range(ga.generationsToConvergence()):
        print("Generation {} is done.".format(ga.atGeneration()))
        ga.advanceToNextGeneration()
        text_file.write(ga.report())
    text_file.close()
    print(ga.report())

    # compute the model responses by temperature for plotting

    P = ga.bestCreatureParam()
    print("The optimal model parameters determined are:")
    print(P)
    model1 = np.zeros(stress1.size)
    for i in range(stress1.size):
        model1[i] = (ss.theta(P, temp[0]) * ss.Zener(P, stress1[i])
                     / math.sqrt(3.0))
    model2 = np.zeros(stress2.size)
    for i in range(stress2.size):
        model2[i] = (ss.theta(P, temp[1]) * ss.Zener(P, stress2[i])
                     / math.sqrt(3.0))
    model3 = np.zeros(stress3.size)
    for i in range(stress3.size):
        model3[i] = (ss.theta(P, temp[2]) * ss.Zener(P, stress3[i])
                     / math.sqrt(3.0))
    model4 = np.zeros(stress4.size)
    for i in range(stress4.size):
        model4[i] = (ss.theta(P, temp[3]) * ss.Zener(P, stress4[i])
                     / math.sqrt(3.0))

    # compute the identity line
    ideal = np.array([1.0e-7, 1.0e-4])

    # create the plot object
    ax = plt.subplot(1, 1, 1)
    ax.set_xscale('log')
    ax.set_yscale('log')

    # select the font family and allow TeX commands
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)

    # change font size for minor and major tick labels
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tick_params(axis='both', which='minor', labelsize=10)

    # add the curves
    line1, = ax.plot(creepRate1, model1, 'yv')
    line2, = ax.plot(creepRate2, model2, 'b^')
    line3, = ax.plot(creepRate3, model3, 'g<')
    line4, = ax.plot(creepRate4, model4, 'r>')
    line5, = ax.plot(ideal, ideal, 'k-')

    # create the plot
    plt.title('Steady-State Creep', fontsize=20)
    plt.xlabel('Experimental Creep Rate, 1/s', fontsize=16)
    plt.ylabel('Modeled Creep Rate, 1/s', fontsize=16)
    plt.legend([line1, line2, line3, line4, line5],
               ['1089 K', '1255 K', '1422 K', '1589 K', 'agreement'],
               bbox_to_anchor=(0.05, 0.925), loc=2, borderaxespad=0.0)
    plt.savefig('ssCreep')
    plt.show()


runGenAlg()
