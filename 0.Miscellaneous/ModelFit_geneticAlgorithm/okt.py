from numpy                import array, sum, inf, square, mean, sqrt, var, logspace, log10
from tkinter              import Tk
from matplotlib.pyplot    import figure, tight_layout, legend, yscale, xscale, savefig, scatter, plot, title
from matplotlib           import rcParams, rcParamsDefault
from pandas               import read_csv
from tkinter              import filedialog
from tkinter.simpledialog import askinteger
from warnings             import filterwarnings
from sys                  import exit as sysexit
from scipy                import optimize
from scipy.optimize       import differential_evolution

def endProgram():
    sysexit()

def setupTkinterWindow():
    root = Tk()
    root.withdraw()

def askForCVSfile():
    filename = filedialog.askopenfilename(filetypes=[("CSV files", ".csv")])
    return filename

def askForNumberOfMaxwellElements():
    numberOfMaxwellElements = askinteger(
        'Enter the number of Maxwell elements',
        'Enter the number of Maxwell elements.\nEntering a negative value OR canceling will end the program.')
    return numberOfMaxwellElements

def getDFfromCSV():
    filename = askForCVSfile()
    df = read_csv(filename, delimiter=",")
    return df

def readCSVfile():
    try:
        return getDFfromCSV()
    except:
        input("File search was canceled OR the selected file is not a CSV.\nPress ENTER to close.")
        endProgram()

def initializePlotAndSetPlotSize():
    scale  = 6
    plotFigure = figure(figsize=(3*scale, 2*scale))
    tight_layout()
    return plotFigure

def getCurrentAxesInstance(plotFigure):
    return plotFigure.gca()

def showPlotLegend(ax):
    handles, labels = ax.get_legend_handles_labels();
    lgd = dict(zip(labels, handles))
    legend(lgd.values(), lgd.keys(), prop={'size': 22}, loc="best")
    
def setPlotTitle(numberOfMaxwellElements):
    title(str(numberOfMaxwellElements) + " Maxwell Elements", size=24);

def namePlotAxes(ax):
    ax.set_xlabel("$\omega$", fontsize=24)
    ax.set_ylabel("${G^{\prime\prime}(\omega)}$", fontsize=24)

def formatTicksAndLabelFontSizes(ax):
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(18)
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(18)
    ax.tick_params(which='both', direction='in', length=5, width=2, bottom=True, top=True, left=True, right=True)
    
def setLogScale():
    yscale('log')
    xscale('log')
    
def savePlot(numberOfMaxwellElements):
    savefig(
        str(numberOfMaxwellElements) + 'elements_DiscreteRelaxationSpectra.png',
        dpi=200,
        bbox_inches='tight')
    
def recoverMatplotlibDefaults():
    rcParams.update(rcParamsDefault)

def getColumnNames(df):
    x_str = df.columns[0]
    y_str = df.columns[1]
    return x_str, y_str

def removeNANsFromData(df, x_str, y_str):
    df = df.dropna(subset=[x_str, y_str])
    return df
    
def getXYdata(df):
    x_str, y_str = getColumnNames(df)
    df           = removeNANsFromData(df, x_str, y_str)
    x            = df[x_str];
    y            = df[y_str];
    x, y         = pandasSeries2numpyArray(x, y)
    return x, y

def scatterExperimentalData(x, y):
    scatter(x, y, s=45, marker='o', label="data")
    plot(x, y, linewidth=1, linestyle='-.')

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

def generatePlot(x, y, xModel, yModel, numberOfMaxwellElements):
    plotFigure = initializePlotAndSetPlotSize()
    ax0        = getCurrentAxesInstance(plotFigure)
    scatterExperimentalData(x, y)
    plotFittedCurve(xModel, yModel)
    namePlotAxes(ax0)
    formatTicksAndLabelFontSizes(ax0)
    setLogScale()
    showPlotLegend(ax0)
    setPlotTitle(numberOfMaxwellElements)
    savePlot(numberOfMaxwellElements)
    recoverMatplotlibDefaults()

def evaluateNumberOfMaxwellElements(numberOfMaxwellElements):
    if str(type(numberOfMaxwellElements)) == "<class 'NoneType'>":
        endProgram()
    elif numberOfMaxwellElements < 1:
        print("The number of Maxwell elements shall be greater than zero.")
        endProgram()
        
def pandasSeries2numpyArray(x, y):
    x = array(x)
    y = array(y)
    return x, y

def getMaxValue(x, y):
    maxX  = max(x)
    maxY  = max(y)
    maxXY = max(maxX, maxY)
    return maxXY

def generate_Initial_Parameters(x, y, numberOfMaxwellElements):
    # min and max used for bounds
    maxX  = max(x)
    maxY  = max(y)
    maxXY = max(maxX, maxY)
    
    parameterBounds = []
    for i in range(numberOfMaxwellElements*2):
        parameterBounds.append([0, maxXY])
    
    # function for genetic algorithm to minimize (sum of squared error)
    def sumOfSquaredError(parameterTuple):
        filterwarnings("ignore") # do not print warnings by genetic algorithm
        val = Maxwell_lossModuli(x, *parameterTuple)
        return sum((y - val) ** 2.0)
    # "seed" the numpy random number generator for repeatable results
    result = differential_evolution(sumOfSquaredError, parameterBounds, seed=3)
    return result.x
    
def fitdata(x, y):
    numberOfMaxwellElements = askForNumberOfMaxwellElements()
    evaluateNumberOfMaxwellElements(numberOfMaxwellElements)
    geneticParameters = generate_Initial_Parameters(x, y, numberOfMaxwellElements)

    # curve fit the test data
    numberOfParameters = len(geneticParameters)
    upbound = [inf]*numberOfParameters
    
    try:
        fittedParameters, pcov = optimize.curve_fit(Maxwell_lossModuli, x, y, geneticParameters, bounds=(0, upbound));
    except:
        print("Optimal parameters not found: The maximum number of function evaluations is exceeded.\nTry a different number of Maxwell elements (or add more data points).")
        return 0, 0, numberOfMaxwellElements

    modelPredictions = Maxwell_lossModuli(x, *fittedParameters) 

    absError = modelPredictions - y

    SE = square(absError) # squared errors
    MSE = mean(SE) # mean squared errors
    RMSE = sqrt(MSE) # Root Mean Squared Error, RMSE
    Rsquared = 1.0 - (var(absError) / var(y))
    
    if Rsquared < 0.85:
        print("Frequency range is not enough to specify " + str(numberOfMaxwellElements) + " relaxation points \nTry a different number of Maxwell elements (or add more data points).")
        return 0, 0, numberOfMaxwellElements
    
    print(">>> " + str(numberOfMaxwellElements) + " Maxwell Elements")
    
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

    minxMagnitude = log10(min(x))
    maxxMagnitude = log10(max(x))
    xModel = logspace(minxMagnitude, maxxMagnitude)
    yModel = Maxwell_lossModuli(xModel, *fittedParameters)
    
    return xModel, yModel, numberOfMaxwellElements

def plotFittedCurve(xModel, yModel):
    plot(xModel, yModel, linestyle='-', linewidth=3, label="fit")