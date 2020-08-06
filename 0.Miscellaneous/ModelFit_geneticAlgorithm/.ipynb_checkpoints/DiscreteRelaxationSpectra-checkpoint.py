%matplotlib inline
import okt
import sys

df             = okt.readCSVfile()
sys.exit()
x, y           = okt.getXYdata(df)
figure         = okt.initializePlotAndSetPlotSize()
ax0            = okt.getCurrentAxesInstance(figure)
xModel, yModel = okt.fitdata(x, y)
okt.scatterExperimentalData(x, y)
okt.plotFittedCurve(xModel, yModel)
okt.namePlotAxes(ax0)
okt.formatTicksAndLabelFontSizes(ax0)
okt.setLogScale()
okt.showPlotLegend(ax0)
okt.saveAndShowPlot()
okt.recoverMatplotlibDefaults()