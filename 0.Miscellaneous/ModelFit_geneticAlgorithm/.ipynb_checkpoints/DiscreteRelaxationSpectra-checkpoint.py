import okt

okt.setupTkinterWindow()
df   = okt.readCSVfile()
x, y = okt.getXYdata(df)
    
while(1):
    model = okt.fitdata(x, y)
    okt.generatePlot(x, y, model[0], model[1], model[2])
    print()