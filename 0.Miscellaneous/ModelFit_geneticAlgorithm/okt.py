from pandas import read_csv
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl
from warnings import filterwarnings
from sys import exit as sysexit

from scipy import optimize
from scipy.optimize import differential_evolution

def askForCSVfile(filename):
    df = read_csv(filename, delimiter=",")
    return df

def readCSVfile(filename):
    try:
        return askForCSVfile(filename)
    except:
        sysexit("File search was canceled OR Selected file is not a CSV")

def initializePlotAndSetPlotSize():
    scale  = 6
    figure = plt.figure(figsize=(3*scale, 2*scale))
    plt.tight_layout()
    return figure

def getCurrentAxesInstance(figure):
    return figure.gca()

def showPlotLegend(ax):
    handles, labels = ax.get_legend_handles_labels();
    lgd = dict(zip(labels, handles))
    plt.legend(lgd.values(), lgd.keys(), prop={'size': 22}, loc="best")
    
def namePlotAxes(ax):
    ax.set_xlabel("$\omega$", fontsize=24)
    ax.set_ylabel("${G^{\prime\prime}(\omega)}$", fontsize=24)

def formatTicksAndLabelFontSizes(ax):
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(18)
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(18)
    ax.tick_params(which='both', direction='in', length=5, width=2, bottom=True, top=True, left=True, right=True)
    
def setLogScale():
    plt.yscale('log')
    plt.xscale('log')
    
def saveAndShowPlot():
    plt.savefig('DiscreteRelaxationSpectra.png', dpi=200, bbox_inches='tight')
    plt.show()
    
def recoverMatplotlibDefaults():
    mpl.rcParams.update(mpl.rcParamsDefault)

def getColumnNames(df):
    x_str = df.columns[0]
    y_str = df.columns[1]
    return x_str, y_str

def removeNANsFromData(df, x_str, y_str):
    df = df.dropna(subset=[x_str, y_str])
    return df
    
def getXYdata(df):
    x_str, y_str = getColumnNames(df)
    df = removeNANsFromData(df, x_str, y_str)
    x = df[x_str];
    y = df[y_str];
    return x, y

def scatterExperimentalData(x, y):
    plt.scatter(x, y, s=45, marker='o', label="data")
    plt.plot(x, y, linewidth=1, linestyle='-.')

def Maxwell_lossModuli(omega_, *p):
    eta_    = p[0            :int(len(p)/2)]
    lambda_ = p[int(len(p)/2):len(p)       ]
    
    sum_ = 0
    for i in range(len(eta_)):
        nume = eta_[i] * omega_
        deno = 1 + (omega_ * lambda_[i])**2
        res  = nume/deno
        sum_ = sum_ + res
    
    return sum_

def fitdata(x, y, numberOfMaxwellElements):
    if numberOfMaxwellElements < 1:
        sysexit("The number of Maxwell elements shall be greater than zero.")
    
    x = np.array(x)
    y = np.array(y)
    
    # function for genetic algorithm to minimize (sum of squared error)
    def sumOfSquaredError(parameterTuple):
        filterwarnings("ignore") # do not print warnings by genetic algorithm
        val = Maxwell_lossModuli(x, *parameterTuple)
        return np.sum((y - val) ** 2.0)
    
    def generate_Initial_Parameters(numberOfMaxwellElements):
        # min and max used for bounds
        maxX = max(x)
        maxY = max(y)
        maxXY = max(maxX, maxY)

        parameterBounds = []
        
        for i in range(numberOfMaxwellElements*2):
            parameterBounds.append([0, maxXY])
        
        # "seed" the numpy random number generator for repeatable results
        result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
        return result.x

    # generate initial parameter values
    geneticParameters = generate_Initial_Parameters(numberOfMaxwellElements)

    #print('geneticParameters', geneticParameters)
    
    # curve fit the test data
    numberOfParameters = len(geneticParameters)
    upbound = [np.inf]*numberOfParameters
    
    try:
        fittedParameters, pcov = optimize.curve_fit(Maxwell_lossModuli, x, y, geneticParameters, bounds=(0, upbound));
    except:
        sysexit("Optimal parameters not found: The maximum number of function evaluations is exceeded.\nTry a different number of Maxwell elements")

    modelPredictions = Maxwell_lossModuli(x, *fittedParameters) 

    absError = modelPredictions - y

    SE = np.square(absError) # squared errors
    MSE = np.mean(SE) # mean squared errors
    RMSE = np.sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (np.var(absError) / np.var(y))
    
    if Rsquared < 0:
        sysexit("Range is not enough to specify " + str(int(numberOfParameters/2)) + " relaxation points \nTry a different number of Maxwell elements (or add more data points)")
    
    print("        n                       ")
    print("      =====                     ")
    print("      \         eta  . omega    ")
    print("       \           i            ")
    print("G'' =   \   --------------------")
    print("        /            2         2")
    print("       /    1 + omega  . lambda ")
    print("      /                        i")
    print("      =====                     ")
    print("      i = 1                     ")
    
    print('eta_i:                \n', fittedParameters[0                        :int(numberOfParameters/2)])
    print('lambda_i:             \n', fittedParameters[int(numberOfParameters/2):numberOfParameters       ])
    print('Root Mean Squared Error:', RMSE)
    print('R-squared:              ', Rsquared)

    #print()
    
    #xModel = np.logspace(min(x), max(x), 100)
    xModel = np.logspace(-2, 3)
    yModel = Maxwell_lossModuli(xModel, *fittedParameters)
    
    return xModel, yModel

def plotFittedCurve(xModel, yModel):
    plt.plot(xModel, yModel, linestyle='-', linewidth=3, label="fit")