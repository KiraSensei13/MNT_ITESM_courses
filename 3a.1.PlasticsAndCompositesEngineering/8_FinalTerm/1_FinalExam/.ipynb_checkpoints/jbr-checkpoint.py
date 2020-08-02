import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.linear_model import LinearRegression
from scipy import special, optimize

#####   #####   #####   #####   #####

# Format scientific notation
def format_e(n):
    a = '%E' % n
    numeric_part    = a.split('E')[0].rstrip('0').rstrip('.')
    int_part        = numeric_part.split('.')[0]
    decimal_part    = numeric_part.split('.')[1]
    scientific_part = a.split('E')[1]
    return int_part + '.' + decimal_part[0:2] + 'e' + scientific_part

#####   #####   #####   #####   #####

def tand(storageModulus=1, lossModulus=1):
    tand = lossModulus / storageModulus;
    return tand;

def calculate_tand(dataframe=pd.DataFrame()):
    sample_id = '';
    for i in range(len(dataframe.columns)):
        if ("$\omega" in dataframe.columns[i]):
            sample_id = dataframe.columns[i].split('_')[1];

            # Compute tan d
            dataframe["${tan\delta}_" + sample_id] = \
            tand(
                dataframe["${G^{\prime}(\omega)}_" + sample_id],
                dataframe["${G^{\prime\prime}(\omega)}_" +sample_id]);

    dataframe = dataframe.reindex(sorted(dataframe.columns), axis=1);
    return dataframe

#####   #####   #####   #####   #####

def complexModulus(storageModulus=1, lossModulus=1):
    Gc = np.sqrt(lossModulus**2 + storageModulus**2);
    return Gc;

def calculate_complexModulus(dataframe=pd.DataFrame()):
    sample_id = '';
    for i in range(len(dataframe.columns)):
        if ("$\omega" in dataframe.columns[i]):
            sample_id = dataframe.columns[i].split('_')[1];

            # Compute tan d
            dataframe["${G^{*}(\omega)}_" + sample_id] = \
            complexModulus(
                dataframe["${G^{\prime}(\omega)}_" + sample_id],
                dataframe["${G^{\prime\prime}(\omega)}_" +sample_id]);

    dataframe = dataframe.reindex(sorted(dataframe.columns), axis=1);
    return dataframe

#####   #####   #####   #####   #####

def calculate_shiftedCurves(dataframe=pd.DataFrame(), aT=[], bT=[]):
    sample_id = '';
    temperature_cnt = 0;
    for i in range(len(dataframe.columns)):
        if ("$\omega" in dataframe.columns[i]):
            sample_id = dataframe.columns[i].split('_')[1];

            dataframe["$a_T \omega_" + sample_id] = dataframe["$\omega_" + sample_id] * aT[temperature_cnt]
            dataframe["$b_T {G^{*}(\omega)}_" + sample_id] = dataframe["${G^{*}(\omega)}_" + sample_id] * bT[temperature_cnt]
            dataframe["$b_T {G^{\prime\prime}(\omega)}_" + sample_id] = dataframe["${G^{\prime\prime}(\omega)}_" + sample_id] * bT[temperature_cnt]

            temperature_cnt = temperature_cnt + 1

    dataframe = dataframe.reindex(sorted(dataframe.columns), axis=1);
    return dataframe

#####   #####   #####   #####   #####

def _invT(T, T0):
    invT = (1/T) - (1/T0);
    return invT;

def calculate_activationEnergies(T=[], Tref=1, aT=[], bT=[]):
    # Set plot size and axis labels' font size
    pltname = "slope of the lin. fit is $E_H$ & $E_V$";
    scale   = 6;
    fig     = plt.figure(figsize=(3*scale, 2*scale));
    plt.rc('xtick', labelsize=15);
    plt.rc('ytick', labelsize=15);
    plt.tight_layout();
    ax0 = plt.gca();

    # calculate varibles
    x = _invT(T, Tref);

    # HORIZONTAL ACTIVATION ENERGY EH

    # perform a linear regression
    y = np.log(aT);
    model = LinearRegression().fit(np.array(x).reshape((-1, 1)), np.array(y));
    # get the slope
    E_H = model.coef_[0];

    # plot variables and lin. fit
    plt.scatter(x, y, s=25, label='horizontal activation energy ' + r'$E_H$' + " = " + str(round(E_H, 2)) + r'$\frac{kJ}{mol}$');
    plt.plot(x, model.predict(np.array(x).reshape((-1, 1))), linewidth=1);

    # VERTICAL ACTIVATION ENERGY EV

    # perform a linear regression
    y = np.log(bT);
    model = LinearRegression().fit(np.array(x).reshape((-1, 1)), np.array(y));
    # get the slope
    E_V = model.coef_[0];

    # plot variables and lin. fit
    plt.scatter(x, y, s=25, label='vertical activation energy ' + r'$E_V$' + " = " + str(round(E_V, 2)) + r'$\frac{kJ}{mol}$');
    plt.plot(x, model.predict(np.array(x).reshape((-1, 1))), linewidth=1);

    # Show the plot lengend to link colors and polymer names
    handles, labels = ax0.get_legend_handles_labels();
    lgd = dict(zip(labels, handles));

    # fig.autofmt_xdate();
    ax0.set_xlabel(r'$\frac{1}{T} - \frac{1}{T_0}$', fontsize=24);
    ax0.set_ylabel(r'$Log(a_T)$' + ' and ' + r'$Log(b_T)$', fontsize=24);

    for tick in ax0.xaxis.get_major_ticks(): tick.label.set_fontsize(18);
    for tick in ax0.yaxis.get_major_ticks(): tick.label.set_fontsize(18);
    ax0.tick_params(which='both', direction='in', length=5, width=2, bottom=True, top=True, left=True, right=True);

    # Display main plot
    plt.yscale('linear');
    plt.xscale('linear');
    plt.legend(lgd.values(), lgd.keys(), prop={'size': 18}, loc="best");
    plt.title(pltname, size=24);
    plt.savefig(pltname + '.png', dpi=200, bbox_inches='tight');
    plt.show();
    mpl.rcParams.update(mpl.rcParamsDefault); # Recover matplotlib defaults
    
    return E_H, E_V
    
#####   #####   #####   #####   #####

def calculate_shiftFactorWithWLF(E=1, T=1, Tref=1):
    return np.exp(E*_invT(T, Tref))    

#####   #####   #####   #####   #####

def complexViscosity(complexModulus=1, angularFrequency=1):
    etac = complexModulus/angularFrequency;
    return etac;

def calculate_complexViscosity(dataframe=pd.DataFrame()):
    sample_id = '';
    for i in range(len(dataframe.columns)):
        if ("$\omega" in dataframe.columns[i]):
            sample_id = dataframe.columns[i].split('_')[1];

            # Compute
            dataframe["${\eta^{*}(\omega)}_" + sample_id] = \
            complexViscosity(
                dataframe["$b_T {G^{*}(\omega)}_" + sample_id],
                dataframe["$a_T \omega_" + sample_id]);

    dataframe = dataframe.reindex(sorted(dataframe.columns), axis=1);
    return dataframe

#####   #####   #####   #####   #####

def get_masterCurveData(df=pd.DataFrame(), tempStr=210, a=1, b=1):

    master_omega_G2 = df["$a_T \omega_{170°C}$"] * a
    master_omega_G2 = master_omega_G2.append(df["$a_T \omega_{190°C}$"] * a)
    master_omega_G2 = master_omega_G2.append(df["$a_T \omega_{210°C}$"] * a)
    master_omega_G2 = master_omega_G2.append(df["$a_T \omega_{230°C}$"] * a)
    
    master_omega_eta = master_omega_G2.append(df["$\dot{\gamma}_{210°C}$"] * a)

    master_eta = df["${\eta^{*}(\omega)}_{170°C}$"] * b
    master_eta = master_eta.append(df["${\eta^{*}(\omega)}_{190°C}$"] * b)
    master_eta = master_eta.append(df["${\eta^{*}(\omega)}_{210°C}$"] * b)
    master_eta = master_eta.append(df["${\eta^{*}(\omega)}_{230°C}$"] * b)
    master_eta = master_eta.append(df["${\eta(\dot{\gamma})}_{210°C}$"] * b)
    
    master_G2 = df["$b_T {G^{\prime\prime}(\omega)}_{170°C}$"] * b
    master_G2 = master_G2.append(df["$b_T {G^{\prime\prime}(\omega)}_{190°C}$"] * b)
    master_G2 = master_G2.append(df["$b_T {G^{\prime\prime}(\omega)}_{210°C}$"] * b)
    master_G2 = master_G2.append(df["$b_T {G^{\prime\prime}(\omega)}_{230°C}$"] * b)

    return master_omega_eta, master_eta, master_omega_G2, master_G2

#####   #####   #####   #####   #####

def _alpha(n_1, lambda_, dot_gamma_o):
    nume = 1 + n_1*lambda_*dot_gamma_o
    deno = lambda_
    res = nume / deno
    return res

def _beta(n_2, lambda_, dot_gamma_o):
    nume = 1 + n_2*lambda_*dot_gamma_o
    deno = lambda_
    res = nume / deno
    return res

def Wagner_eta_infty(dot_gamma, *p):
    a_          = p[0:5]
    lambda_     = p[5:10]
    f_1         = p[10]
    f_2         = 1 - f_1
    n_1         = p[11]
    n_2         = p[12]
    
    sum_1 = 0
    for i in range(0, 5, 1):
        alpha = _alpha(n_1, lambda_[i], dot_gamma)
        res_1 = a_[i] / alpha**2
        sum_1 = sum_1 + res_1
        
    sum_2 = 0
    for i in range(0, 5, 1):
        beta  = _beta(n_2, lambda_[i], dot_gamma)
        res_2 = a_[i] / beta**2
        sum_2 = sum_2 + res_2
    
    res = (f_1 * sum_1) + (f_2 * sum_2)
    
    return res/10

def Wagner_fit_eta(t, eta):    
    # Initial guess
    ai      = [35000, 11000, 4000, 1000, 20]
    lambdai = [0.00005, 0.005, 0.05, 0.5, 5]
    f1 = [0.5]
    n1 = [2.8]
    n2 = [0.07]
    p = ai + lambdai + f1 + n1 + n2
    
    # Fit the model
    upbound    = [np.inf]*10 + [1] + [np.inf]*2  
    model      = optimize.curve_fit(Wagner_eta_infty, t, eta, p) #, bounds=(0, upbound)); #bounds=(0, [3., 1., 0.5])
    parameters = model[0]

    # Show the fitting parameters
    print(parameters)
    
    return parameters

#####   #####   #####   #####   #####

def plot_Wagner_fit(df_master=pd.DataFrame(), x_str="", y_str=""):

    parameters = Wagner_fit_eta(
        pd.Series(df_master[x_str]).dropna(),
        pd.Series(df_master[y_str]).dropna())
    
    # Set plot size and axis labels' font size
    pltname = 'Wagner fit';
    scale   = 6;
    fig     = plt.figure(figsize=(3*scale, 2*scale));
    plt.rc('xtick', labelsize=15);
    plt.rc('ytick', labelsize=15);
    plt.tight_layout();
    ax0 = plt.gca();
    
    # Stablish the plot area
    ax0 = plt.gca()

    # plot dataset
    t   = pd.Series(df_master[x_str]).dropna()
    eta = pd.Series(df_master[y_str]).dropna()
    plt.scatter(t, eta, label='experimental')
    plt.plot(t, eta, linewidth=1, linestyle=':')

    # Plot fit
    n_1 = parameters[11]
    n_2 = parameters[12]
    t = np.logspace(-2, 4, 100)
    eta_fit = Wagner_eta_infty(t, *parameters)
    plt.plot(t, eta_fit, linewidth=3, label = r'$f_1 = $' + str(round(parameters[10],2)) + ', ' +
             r'$n_1 = $' + str(round(n_1, 2)) + ', ' +
             r'$n_2 = $' + str(round(n_2, 2)));

    # Format and Display plots
    ax0.tick_params(which='both', direction='in', width=2, bottom=True, top=True, left=True, right=True);
    plt.yscale('log');
    plt.xscale('log');

    # Show the plot lengend to link colors and polymer names
    handles, labels = ax0.get_legend_handles_labels();
    lgd = dict(zip(labels, handles));
    
    plt.xlabel(r'$1/s$', fontsize=24);
    plt.ylabel(r'$Pa \cdot s$', fontsize=24);
    plt.legend(lgd.values(), lgd.keys(), prop={'size': 22}, loc="best");
    plt.title(pltname, size=24);
    plt.savefig(pltname + '.png', dpi=200, bbox_inches='tight');
    plt.show();
    mpl.rcParams.update(mpl.rcParamsDefault); # Recover matplotlib defaults
    
    return parameters

#####   #####   #####   #####   #####

def _N1(dot_gamma, *p):
    a_          = p[0:5]
    lambda_     = p[5:10]
    f_1         = p[10]
    f_2         = 1 - f_1
    n_1         = p[11]
    n_2         = p[12]
    
    sum_1 = 0
    for i in range(0, 5, 1):
        alpha = _alpha(n_1, lambda_[i], dot_gamma)
        res_1 = a_[i] * (alpha**3)
        sum_1 = sum_1 + res_1
        
    sum_2 = 0
    for i in range(0, 5, 1):
        beta  = _beta(n_2, lambda_[i], dot_gamma)
        res_2 = a_[i] * (beta**3)
        sum_2 = sum_2 + res_2
    
    res = dot_gamma**2 * (f_1*sum_1 + f_2*sum_2)
    
    return res/10

def plot_steadyStateN1(parameters):
    # Draw plot canvas
    plotname = 'Steady State 1st Normal Stress Difference';
    scale   = 6;
    fig     = plt.figure(figsize=(3*scale, 2*scale));
    plt.rc('xtick', labelsize=15);
    plt.rc('ytick', labelsize=15);
    plt.tight_layout();
    ax0 = plt.gca();

    a_          = parameters[0:5]
    lambda_     = parameters[5:10]
    f_1         = parameters[10]
    f_2         = 1 - f_1
    n_1         = parameters[11]
    n_2         = parameters[12]

    # Plot fit
    gamma = np.logspace(-2, 3, 100)
    N1 = _N1(gamma, *parameters)
    plt.plot(gamma, N1, linestyle='--', linewidth=3, label = 
             r'$a_1 = $' + format_e(a_[0]) + ", " +
             r'$a_2 = $' + format_e(a_[1]) + ", " +
             r'$a_3 = $' + format_e(a_[2]) + ",\n" +
             r'$a_4 = $' + format_e(a_[3]) + ", " +
             r'$a_5 = $' + format_e(a_[4]) + ",\n" +
             r'$\lambda_1 = $' + format_e(lambda_[0]) + ", " +
             r'$\lambda_2 = $' + format_e(lambda_[1]) + ", " +
             r'$\lambda_3 = $' + format_e(lambda_[2]) + ",\n" +
             r'$\lambda_4 = $' + format_e(lambda_[3]) + ", " +
             r'$\lambda_5 = $' + format_e(lambda_[4]) + ",\n" +
             r'$f_1 = $' + format_e(f_1) + ", " +
             r'$f_2 = $' + format_e(f_2) + ",\n" +
             r'$n_1 = $' + format_e(n_1) + ", " +
             r'$n_2 = $' + format_e(n_2));
    
    gammas = [0.1, 1, 10, 100]
    for g in gammas:
        N1 = _N1(g, *parameters)
        plt.scatter(g, N1, s=90, label=r'$N1(' + str(g) + ') = ' + format_e(N1) + '$')

    # Format and Display plots
    ax0.tick_params(which='both', direction='in', width=2, bottom=True, top=True, left=True, right=True);
    plt.yscale('log');
    plt.xscale('log');
    plt.xlabel(r'$1/s$', fontsize=24);
    plt.ylabel(r'$Pa$', fontsize=24);
    plt.title(plotname, size=24);
    plt.legend(prop={'size': 22});
    plt.savefig('plt_' + plotname + '.png', dpi=200, bbox_inches='tight');
    plt.show();
    mpl.rcParams.update(mpl.rcParamsDefault); # Recover matplotlib defaults

#####   #####   #####   #####   #####

def plot(dataframe=pd.DataFrame(), pltname="plotTitle", x_str=[], x_units='', y_str=[], y_units=''):
    # Set plot size and axis labels' font size
    scale = 6;
    fig   = plt.figure(figsize=(3*scale, 2*scale));
    plt.rc('xtick', labelsize=15);
    plt.rc('ytick', labelsize=15);
    plt.tight_layout();

    # Stablish the plot area
    ax0 = plt.gca();

    for xName, yName in zip(x_str, y_str):
        # Remove NANs from interesting x,y data
        df_fil = pd.DataFrame(dataframe);
        df_fil = df_fil.dropna(subset=[xName, yName]);

        # Extract data from a specific country
        x = df_fil.iloc[:][xName];
        y = df_fil.iloc[:][yName];

        # Scatter the data and plot a curve to join the points
        plt.scatter(x, y, s=45, marker='o', label=yName);
        plt.plot(x, y, linewidth=1, linestyle='-.');

    # Show the plot lengend to link colors and polymer names
    handles, labels = ax0.get_legend_handles_labels();
    lgd = dict(zip(labels, handles));

    # fig.autofmt_xdate();
    ax0.set_xlabel(str(x_units), fontsize=24);
    ax0.set_ylabel(str(y_units), fontsize=24);

    for tick in ax0.xaxis.get_major_ticks(): tick.label.set_fontsize(18);
    for tick in ax0.yaxis.get_major_ticks(): tick.label.set_fontsize(18);
    ax0.tick_params(which='both', direction='in', length=5, width=2, bottom=True, top=True, left=True, right=True);

    # Display main plot
    plt.yscale('log');
    plt.xscale('log');
    plt.legend(lgd.values(), lgd.keys(), prop={'size': 22}, loc="best");
    plt.title(pltname, size=24);
    plt.savefig(pltname + '.png', dpi=200, bbox_inches='tight');
    plt.show();
    mpl.rcParams.update(mpl.rcParamsDefault); # Recover matplotlib defaults