import matplotlib.pyplot as plt
import numpy as np
import math
import random


class reaction:
    # Print out help statement in case of incorrect definition or syntax
    def helpPrint():
        print("Please specify a reaction in the form:")
        print("\n  Reactants -> Products")
        print("\nSome possible reactions include:")
        print("  - 1st order (monomolecular): A -> B")
        print("  - 2nd order (bimolecular):  2A -> B,  A + B -> C")
        print("  - 3rd order (trimolecular): 3A -> B,  A + 2B -> C,  A + B + C -> D")

    # Construct chemical reaction class instance
    def __init__(self, reactionStr, reactionConst):
        self.reactants = []
        self.products = []
        self.coeffsReact = []
        self.coeffsProd = []
        self.rate = reactionConst
        self.reactionStr = ""

        # Check if reaction rate depends on variable

        # If not a valid format for a reaction
        if '->' not in reactionStr:
            print("ERROR: Invalid reaction syntax (\'->\' not in reaction string)")
            helpPrint()
            exit(1)
        else:
            self.reactionStr = reactionStr

            reactProd = reactionStr.split("->")
            if len(reactProd) > 2:
                print("ERROR: Invalid reaction syntax (\'->\' appears more than once in reaction string)")
                helpPrint()
                exit(1)
            if reactProd[0].strip() == "":
                print("ERROR: Reaction has no reactants")
                helpPrint()
                exit(1)
            else:
                reactStrs = reactProd[0].split(' + ')
                for reactStr in reactStrs:
                    currCoeffStr = ""
                    currMolecStr = ""
                    charsFound = False
                    for char in reactStr:
                        if char.isdigit() and charsFound == False:
                            currCoeffStr += char
                        elif char.isspace():
                            continue
                        else:
                            charsFound = True
                            currMolecStr += char
                    self.reactants.append(currMolecStr)
                    if currCoeffStr == "":
                        self.coeffsReact.append(1)
                    else:
                        self.coeffsReact.append(int(currCoeffStr))

            if reactProd[1].strip() == "":
                print("ERROR: Reaction has no products")
                helpPrint()
                exit(1)
            else:
                prodStrs = reactProd[1].split(' + ')
                for prodStr in prodStrs:
                    currCoeffStr = ""
                    currMolecStr = ""
                    charsFound = False
                    for char in prodStr:
                        if char.isdigit() and charsFound == False:
                            currCoeffStr += char
                        elif char.isspace():
                            continue
                        else:
                            charsFound = True
                            currMolecStr += char
                    self.products.append(currMolecStr)
                    if currCoeffStr == "":
                        self.coeffsProd.append(1)
                    else:
                        self.coeffsProd.append(int(currCoeffStr))
    
    # Return string defining the reaction
    def getReactionStr(self):
        return self.reactionStr

    # Return key-value pairs of reactants and their corresponding coefficients
    def getReactants(self):
        reactDict = {}
        for i in range(len(self.reactants)):
            reactDict[self.reactants[i]] = self.coeffsReact[i]
        return reactDict

    # Return key-value pairs of products and their corresponding coefficients
    def getProducts(self):
        prodDict = {}
        for i in range(len(self.products)):
            prodDict[self.products[i]] = self.coeffsProd[i]
        return prodDict

    # Return the propensity for the reaction, given the number of each reactant
    def getPropensity(self, numReactants):
        propensityVect = []
        if len(numReactants) != len(self.reactants):
            print("ERROR: Mismatch between ammount vector and number of reactants")
            return -1.0
        else:
            for i in range(len(numReactants)):
                if numReactants[i] == -1:
                    currProp = 1
                else:
                    currProp = math.comb(numReactants[i], self.coeffsReact[i])
                propensityVect.append(currProp)
            return self.rate * np.prod(propensityVect)




# Implement Gillespie algorithm for discrete, stochastic
# simulation of chemical reactions
#
#   Given:
#     - Set of N molecules and each's discrete number
#     - Set of M reactions and each's reaction constant
#
#      R : List of strings defining reactions in the form '2A -> B' or 'A + B -> C' etc.
#      X : Dictionary of (molecule, count) pairs representing the number of each of the N molecules
#      C : List of floating points representing the reaction constants for M reactions
#   maxT : The amount of time through which to propagate the stochastic simulation
#
def gillespie(R, X, C, maxT, t0 = 0.0):

    if len(R) != len(C):
        print("ERROR: Number of reactions should match number of reaction rates")
        return 0

    # Obtain number of unique molecules
    N = len(X)

    # Obtain number of possible reactions
    M = len(C)

    # Initialize empty time and molecular count vectors
    t_vect = []
    count_dict = {}
    for Xi in X:
        count_dict[Xi] = [X[Xi]]

    # Initialize time at t = 0, number of reactions at n = 0
    t = t0
    t_vect.append(t)

    # Instantiate reactions and store in list
    reactionList = []
    for i in range(len(R)):
        currReaction = reaction(R[i].getReactionStr(), C[i])
        reactionList.append(currReaction)

    # Perform main Gillespie loop
    while t < maxT:
        propVect = []
        for rxn in reactionList:
            currReactCountSet = []
            currReacts = rxn.getReactants()
            for reactant in currReacts:
                currReactCountSet.append(X[reactant])
            currProp = rxn.getPropensity(currReactCountSet)
            propVect.append(currProp)
        totalPropensity = sum(propVect)

        # Generate two random numbers for stochastic simulation
        r1 = random.uniform(0, 1)
        r2 = random.uniform(0, 1)

        # Use r1 to determine wait time before next reaction
        tau = (1 / totalPropensity) * np.log(1 / r1)

        # Use r2 to determine which reaction will occur
        mu = 0
        N = r2 * totalPropensity - propVect[mu]
        while N > 0:
            mu += 1
            N -= propVect[mu]
        nextReaction = reactionList[mu]

        # Step forward in time until the next reaction occurs
        t += tau
        print(t)
        # Alter count of molecule species post-reaction
        for Xi in X:
            if Xi in nextReaction.getReactants() and X[Xi] != -1:
                # If amount of a reactant is -1, it will be held constant
                if X[Xi] != -1:
                    X[Xi] -= nextReaction.getReactants()[Xi]
                else:
                    pass
            if Xi in nextReaction.getProducts():
                X[Xi] += nextReaction.getProducts()[Xi]
            else:
                pass

        # Append time, molecular counts to vectors
        t_vect.append(t)
        for Xi in X:
            count_dict[Xi].append(X[Xi])

    return t_vect, count_dict


def plotReaction(t_vect, count_dict):
    for molecule in count_dict:
        plt.plot(t_vect, count_dict[molecule], label = "# molecules of " + molecule)
    plt.title("Graph of ammounts of reactant & product over time")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # Simulate simple first-order system
    """
    c = 0.5

    react1 = reaction('X -> Z', c)

    counts = {'X': 1000, 'Z': 0}

    t_vect, counts_dict = gillespie([react1], counts, [c], maxT = 10.0)

    plt.figure(1)
    plt.plot(t_vect, counts_dict['X'], label = "# of X")
    plt.plot(t_vect, counts_dict['Z'], label = "# of Z")
    plt.legend()
    plt.title("First order system (X -> Z), with c = 0.5")
    plt.show()
    
    # Simulate two-reaction system where X is held constant
    c1 = 5.0
    c2 = 0.005

    react1 = reaction('X + Y -> 2Y', c1)
    react2 = reaction('2Y -> Z', c2)
    counts = {'X': -1, 'Y': 10, 'Z': 0}

    t_vect, counts_dict = gillespie([react1, react2], counts, [c1, c2], maxT = 10.0)
    
    plt.figure(2)
    plt.plot(t_vect, counts_dict['Y'], label = "# of Y w/ Y0 = 10")
    counts = {'X': -1, 'Y': 3000, 'Z': 0}

    t_vect, counts_dict = gillespie([react1, react2], counts, [c1, c2], maxT = 10.0)
    plt.plot(t_vect, counts_dict['Y'], label = "# of Y w/ Y0 = 3000")

    plt.title("Steady state system (X + Y -> 2Y | 2Y -> Z)")
    plt.legend()
    plt.show()

    # Simulate the Lotka Reaction chemical system
    
    c1 = 10.0
    c2 = 0.01
    c3 = 10.0

    react1 = reaction('X + Y1 -> 2Y1', c1)
    react2 = reaction('Y1 + Y2 -> 2Y2', c2)
    react3 = reaction('Y2 -> Z', c3)

    counts = {'X': -1, 'Y1': 1000, 'Y2': 1000, 'Z': 0}

    t_vect, counts_dict = gillespie([react1, react2, react3], counts, [c1, c2, c3], maxT = 30.0)

    plt.plot(t_vect, counts_dict['Y1'], label = "# of Y1")
    plt.plot(t_vect, counts_dict['Y2'], label = "# of Y2")
    plt.legend()
    plt.title("Lotka chemical system")
    plt.show()
    
    # Simulate the Brusselator system

    #c1 = 5000
    #c2 = 50
    #c3 = 0.00005
    #c4 = 5

    counts = {'X1': -1, 'X2': -1, 'Y1': 1000, 'Y2': 2000, 'Z1': 0, 'Z2': 0}
    a = 0.2

    c1 = 5000
    c2 = 2 * c1 / (a * counts['Y1'])
    c3 = 4 * c1 / (a * counts['Y1']**2 * counts['Y2'])
    c4 = c1 / counts['Y1']

    react1 = reaction('X1 -> Y1', c1)
    react2 = reaction('X2 + Y1 -> Y2 + Z1', c2)
    react3 = reaction('2Y1 + Y2 -> 3Y1', c3)
    react4 = reaction('Y1 -> Z2', c4)

    t_vect, counts_dict = gillespie([react1, react2, react3, react4], counts, [c1, c2, c3, c4], maxT = 10.0)

    plt.figure(1)
    plt.plot(t_vect, counts_dict['Y1'])
    plt.title("Brusselator chemical system (Y1)")
    plt.show()

    plt.figure(2)
    plt.plot(t_vect, counts_dict['Y2'])
    plt.title("Brusselator chemical system (Y2)")
    plt.show()
    
    # Simulate the Oregonator system
    
    counts = {'X1': -1, 'X2': -1, 'X3': -1, 'Y1': 500, 'Y2': 1000, 'Y3': 2000, 'Z1': 0, 'Z2': 0}
    p1 = 2000
    p2 = 50000

    c1 = p1 / counts['Y2']
    c2 = p2 / (counts['Y1'] * counts['Y2'])
    c3 = (p1 + p2) / counts['Y1']
    c4 = 2 * p1 / counts['Y1']**2
    c5 = (p1 + p2) / counts['Y3']
    C = [c1, c2, c3, c4, c5]

    react1 = reaction('X1 + Y2 -> Y1', c1)
    react2 = reaction('Y1 + Y2 -> Z1', c2)
    react3 = reaction('X2 + Y1 -> 2Y1 + Y3', c3)
    react4 = reaction('2Y1 -> Z2', c4)
    react5 = reaction('X3 + Y3 -> Y2', c5)
    R = [react1, react2, react3, react4, react5]

    t_vect, counts_dict = gillespie(R, counts, C, t0 = 1.0, maxT = 3.0)

    plt.figure(1)
    plt.plot(t_vect, counts_dict['Y1'], label = "# of Y1")
    plt.plot(t_vect, counts_dict['Y2'], label = "# of Y2")
    plt.plot(t_vect, counts_dict['Y3'], label = "# of Y3")
    plt.title("Oregonator chemical system")
    plt.legend()
    plt.show()
    """

    # Simulate poleward flux of microtubule
    # Notation: T_e_s
    #  - e (+/-): Polarity of the end
    #  - s (+/-/0): State of the end

    # Define rate constants for reactions
    k1_f = 0.04
    k1_r = 0.1

    k2_f = 0.065
    k2_r = 0.85

    k3_f = 0.009
    k3_r = 0.1

    k4_f = 0.001 * k3_f
    k4_r = k3_r

    k5_f = 3.3
    k5_r = 0.00235

    k6_f = 0.001 * k5_f
    k6_r = k5_r

    k7_f = 0.9
    k7_r = 0.001


    # Reaction 1: Binding/unbinding of minispindles
    react1_for_g = reaction('T_+_+ + Msps -> T_+_+—Msps', k1_f)
    react1_rev_g = reaction('T_+_+—Msps -> T_+_+ + Msps', k1_r)
    react1_for_s = reaction('T_+_- + Msps -> T_+_-—Msps', k1_f)
    react1_rev_s = reaction('T_+_-—Msps -> T_+_- + Msps', k1_r)
    react1_for_n = reaction('T_+_0 + Msps -> T_+_0—Msps', k1_f)
    react1_rev_n = reaction('T_+_0—Msps -> T_+_0 + Msps', k1_r)
    react1 = [react1_for_g, react1_rev_g, react1_for_s, react1_rev_s, react1_for_n, react1_rev_n]

    # Reaction 2: Binding/unbinding of K67A protein
    react2_for_g = reaction('T_+_+—X + K67A -> T_+_+—X—K67A', k2_f)
    react2_rev_g = reaction('T_+_+—X—K67A -> T_+_+—X + K67A', k2_r)
    react2_for_s = reaction('T_+_-—X + K67A -> T_+_-—X—K67A', k2_f)
    react2_rev_s = reaction('T_+_-—X—K67A -> T_+_-—X + K67A', k2_r)
    react2_for_n = reaction('T_+_0—X + K67A -> T_+_0—X—K67A', k2_f)
    react2_rev_n = reaction('T_+_0—X—K67A -> T_+_0—X + K67A', k2_r)
    react2 = [react2_for_g, react2_rev_g, react2_for_s, react2_rev_s, react2_for_n, react2_rev_n]

    # Reaction 3: Binding/unbinding of K59C protein from Mast complex
    react3_for_g = reaction('T_+_+—X[Mast] + K59C -> T_+_+—X[Mast]—K59C', k3_f)
    react3_rev_g = reaction('T_+_+—X[Mast]—K59C -> T_+_+—X[Mast] + K59C', k3_r)
    react3_for_s = reaction('T_+_-—X[Mast] + K59C -> T_+_-—X[Mast]—K59C', k3_f)
    react3_rev_s = reaction('T_+_-—X[Mast]—K59C -> T_+_-—X[Mast] + K59C', k3_r)
    react3_for_n = reaction('T_+_0—X[Mast] + K59C -> T_+_0—X[Mast]—K59C', k3_f)
    react3_rev_n = reaction('T_+_0—X[Mast]—K59C -> T_+_0—X[Mast] + K59C', k3_r)
    react3 = [react3_for_g, react3_rev_g, react3_for_s, react3_rev_s, react3_for_n, react3_rev_n]

    # Reaction 4: Binding/unbinding of K59C protein from general complex
    react4_for_g = reaction('T_+_+—X—Mast + K59C -> T_+_+—X—Mast—K59C', k4_f)
    react4_rev_g = reaction('T_+_+—X—Mast—K59C -> T_+_+—X—Mast + K59C', k4_r)
    react4_for_s = reaction('T_+_-—X—Mast + K59C -> T_+_-—X—Mast—K59C', k4_f)
    react4_rev_s = reaction('T_+_-—X—Mast—K59C -> T_+_-—X—Mast + K59C', k4_r)
    react4_for_n = reaction('T_+_0—X—Mast + K59C -> T_+_0—X—Mast—K59C', k4_f)
    react4_rev_n = reaction('T_+_0—X—Mast—K59C -> T_+_0—X—Mast + K59C', k4_r)
    react4 = [react4_for_g, react4_rev_g, react4_for_s, react4_rev_s, react4_for_n, react4_rev_n]

    # Reaction 5: Binding/unbinding of Mast from K59C complex
    react5_for_g = reaction('T_+_+—X[K59C] + Mast -> T_+_+—X[K59C]—Mast', k5_f)
    react5_rev_g = reaction('T_+_+—X[K59C]—Mast -> T_+_+—X[K59C] + Mast', k5_r)
    react5_for_s = reaction('T_+_-—X[K59C] + Mast -> T_+_-—X[K59C]—Mast', k5_f)
    react5_rev_s = reaction('T_+_-—X[K59C]—Mast -> T_+_-—X[K59C] + Mast', k5_r)
    react5_for_n = reaction('T_+_0—X[K59C] + Mast -> T_+_0—X[K59C]—Mast', k5_f)
    react5_rev_n = reaction('T_+_0—X[K59C]—Mast -> T_+_0—X[K59C] + Mast', k5_r)
    react5 = [react5_for_g, react5_rev_g, react5_for_s, react5_rev_s, react5_for_n, react5_rev_n]

    # Reaction 6: Binding/unbinding of Mast from general complex
    #react6_for_g = reaction('



