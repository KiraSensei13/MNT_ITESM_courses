import okt

filename = "./df_masterCurve_G2.csv"
numberOfMaxwellElements = 5

df             = okt.readCSVfile(filename)
x, y           = okt.getXYdata(df)
xModel, yModel = okt.fitdata(x, y, numberOfMaxwellElements)