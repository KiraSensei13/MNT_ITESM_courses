import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from sklearn.linear_model import LinearRegression
from scipy import special, optimize, signal
from scipy.signal import lfilter

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
            dataframe["$b_T {G^{\prime}(\omega)}_" + sample_id] = dataframe["${G^{\prime}(\omega)}_" + sample_id] * bT[temperature_cnt]

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
    
    master_G1 = df["$b_T {G^{\prime}(\omega)}_{170°C}$"] * b
    master_G1 = master_G1.append(df["$b_T {G^{\prime}(\omega)}_{190°C}$"] * b)
    master_G1 = master_G1.append(df["$b_T {G^{\prime}(\omega)}_{210°C}$"] * b)
    master_G1 = master_G1.append(df["$b_T {G^{\prime}(\omega)}_{230°C}$"] * b)

    return master_omega_eta, master_eta, master_omega_G2, master_G2, master_G1

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

def Maxwell_storageModuli(omega_, *p):
    eta_    = p[0            :int(len(p)/2)]
    lambda_ = p[int(len(p)/2):len(p)       ]
    
    #print(p)
    #print("eta", eta_)
    #print("lambda", lambda_)
    
    sum_ = 0
    for i in range(len(eta_)):
        nume = eta_[i] * lambda_[i] * omega_**2
        deno = 1 + omega_**2 + lambda_[i]**2
        res  = nume/deno
        sum_ = sum_ + res
    
    return sum_

def Maxwell_lossModuli(omega_, *p):
    eta_    = p[0            :int(len(p)/2)]
    lambda_ = p[int(len(p)/2):len(p)       ]
    
    #print(p)
    #print("eta", eta_)
    #print("lambda", lambda_)
    
    sum_ = 0
    for i in range(len(eta_)):
        nume = eta_[i] * omega_
        deno = 1 + omega_**2 + lambda_[i]**2
        res  = nume/deno
        sum_ = sum_ + res
    
    return sum_

#####   #####   #####   #####   #####

def Maxwell_fit_lossModulus(freq, G2, *guess):    
    # Initial guess
    ai      = guess[0:int((len(guess))/2)]          #[35000, 11000, 4000, 1000, 20]
    lambdai = guess[int((len(guess))/2):len(guess)] #[0.00005, 0.005, 0.05, 0.5, 5]
    #f1 = guess[-3]
    #n1 = guess[-2]
    #n2 = guess[-1]
    p = ai + lambdai # + (f1, n1, n2)
    
    # Fit the model
    #upbound    = [np.inf]*10 + [1] + [np.inf]*2
    model      = optimize.curve_fit(Maxwell_lossModuli, freq, G2, p) #, bounds=(0, upbound)); #bounds=(0, [3., 1., 0.5])
    parameters = model[0]

    # Show the fitting parameters
    print(parameters)
    
    return parameters

def Maxwell_fit_storageModulus(freq, G2, *guess):    
    # Initial guess
    ai      = guess[0:int((len(guess))/2)]          #[35000, 11000, 4000, 1000, 20]
    lambdai = guess[int((len(guess))/2):len(guess)] #[0.00005, 0.005, 0.05, 0.5, 5]
    #f1 = guess[-3]
    #n1 = guess[-2]
    #n2 = guess[-1]
    p = ai + lambdai # + (f1, n1, n2)
    
    # Fit the model
    #upbound    = [np.inf]*10 + [1] + [np.inf]*2
    model      = optimize.curve_fit(Maxwell_storageModuli, freq, G2, p) #, bounds=(0, upbound)); #bounds=(0, [3., 1., 0.5])
    parameters = model[0]

    # Show the fitting parameters
    print(parameters)
    
    return parameters

def plot_Maxwell_fit_lossModulus(df_master=pd.DataFrame(), x_str="", y_str="", guess=[]):
    
    #parameters = guess;
    
    parameters = Maxwell_fit_lossModulus(
        pd.Series(df_master[x_str]).dropna(),
        pd.Series(df_master[y_str]).dropna(),
        *guess)
    
    # Set plot size and axis labels' font size
    pltname = 'Maxwell fit - loss modulus';
    scale   = 6;
    fig     = plt.figure(figsize=(3*scale, 2*scale));
    plt.rc('xtick', labelsize=15);
    plt.rc('ytick', labelsize=15);
    plt.tight_layout();
    ax0 = plt.gca();
    
    # Stablish the plot area
    ax0 = plt.gca()

    # plot dataset
    t  = pd.Series(df_master[x_str]).dropna()
    G2 = pd.Series(df_master[y_str]).dropna()
    plt.scatter(t, G2, label='experimental')
    plt.plot(t, G2, linewidth=1, linestyle=':')
    
    # Plot fit
    t = np.logspace(-2, 3, 100)
    eta_fit = Maxwell_lossModuli(t, *parameters)
    plt.plot(t, eta_fit, color='green', linewidth=2, linestyle='--', label = "'Relaxation Spectra.exe' fit");
    
    # Format and Display plots
    ax0.tick_params(which='both', direction='in', width=2, bottom=True, top=True, left=True, right=True);
    plt.yscale('log');
    plt.xscale('log');

    # Show the plot lengend to link colors and polymer names
    handles, labels = ax0.get_legend_handles_labels();
    lgd = dict(zip(labels, handles));
    
    plt.xlabel(r'$\gamma$' + '    ' + r'$1/s$', fontsize=24);
    plt.ylabel(r'$G^{\prime\prime}$' + '    ' + r'$Pa$', fontsize=24);
    plt.legend(lgd.values(), lgd.keys(), prop={'size': 22}, loc="best");
    plt.title(pltname, size=24);
    plt.savefig(pltname + '.png', dpi=200, bbox_inches='tight');
    plt.show();
    mpl.rcParams.update(mpl.rcParamsDefault); # Recover matplotlib defaults
    
    return parameters
    
    
def plot_Maxwell_storageModulus(df_master=pd.DataFrame(), x_str="", y_str="", guess=[]):
    
    parameters = Maxwell_fit_storageModulus(
        pd.Series(df_master[x_str]).dropna(),
        pd.Series(df_master[y_str]).dropna(),
        *guess)
    
    # Set plot size and axis labels' font size
    pltname = 'Maxwell fit - storage modulus';
    scale   = 6;
    fig     = plt.figure(figsize=(3*scale, 2*scale));
    plt.rc('xtick', labelsize=15);
    plt.rc('ytick', labelsize=15);
    plt.tight_layout();
    ax0 = plt.gca();
    
    # Stablish the plot area
    ax0 = plt.gca()

    # plot dataset
    t  = pd.Series(df_master[x_str]).dropna()
    G2 = pd.Series(df_master[y_str]).dropna()
    plt.scatter(t, G2, label='experimental')
    plt.plot(t, G2, linewidth=1, linestyle=':')
    
    # Plot fit
    t = np.logspace(-2, 3, 100)
    eta_fit = Maxwell_storageModuli(t, *parameters)
    plt.plot(t, eta_fit, color='green', linewidth=2, linestyle='--', label = "'Relaxation Spectra.exe' fit");
    
    # Format and Display plots
    ax0.tick_params(which='both', direction='in', width=2, bottom=True, top=True, left=True, right=True);
    plt.yscale('log');
    plt.xscale('log');

    # Show the plot lengend to link colors and polymer names
    handles, labels = ax0.get_legend_handles_labels();
    lgd = dict(zip(labels, handles));
    
    plt.xlabel(r'$\gamma$' + '    ' + r'$1/s$', fontsize=24);
    plt.ylabel(r'$G^{\prime}$' + '    ' + r'$Pa$', fontsize=24);
    plt.legend(lgd.values(), lgd.keys(), prop={'size': 22}, loc="best");
    plt.title(pltname, size=24);
    plt.savefig(pltname + '.png', dpi=200, bbox_inches='tight');
    plt.show();
    mpl.rcParams.update(mpl.rcParamsDefault); # Recover matplotlib defaults
    
#####   #####   #####   #####   #####
    
def plot_Wagner_fit_viscosity(df_master=pd.DataFrame(), x_str="", y_str=""):

    parameters = Wagner_fit_eta(
        pd.Series(df_master[x_str]).dropna(),
        pd.Series(df_master[y_str]).dropna())
    
    #parameters = [210.42, 957.75, 4243.12, 11511.83, 31232.21,
    #     0.000526, 0.005263, 0.052632, 0.52631, 5.263158]
    
    # Set plot size and axis labels' font size
    pltname = 'Wagner fit - viscosity';
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
    t = np.logspace(-2, 4, 100)
    eta_fit = Wagner_eta_infty(t, *parameters)
    plt.plot(t, eta_fit, linewidth=3, label = "Scipy fit: " + r'$f_1 = $' + str(round(parameters[10],2)) + ', ' +
             r'$n_1 = $' + str(round(parameters[11], 2)) + ', ' +
             r'$n_2 = $' + str(round(parameters[12], 2)));

    #JBRtoolParameters = parameters.copy()
    #JBRtoolParameters[0:5]  = [ 2.59334433e+08, 1.18280168e+06, 2.48532480e+05, 9.76160818e+03, 1.65786734e+02]
    #JBRtoolParameters[5:10] = [10.35724846e-03, 1.65264506e-01, 1.83612459e-01, 6.21257664e+00, 8.16912971e+01]
    #JBReta_fit = Wagner_eta_infty(t, *JBRtoolParameters)
    #plt.plot(t, JBReta_fit, linewidth=2, linestyle='--', label = "'Relaxation Spectra.exe' fit");
    
    # Format and Display plots
    ax0.tick_params(which='both', direction='in', width=2, bottom=True, top=True, left=True, right=True);
    plt.yscale('log');
    plt.xscale('log');

    # Show the plot lengend to link colors and polymer names
    handles, labels = ax0.get_legend_handles_labels();
    lgd = dict(zip(labels, handles));
    
    plt.xlabel(r'$\gamma$' + '    ' + r'$1/s$', fontsize=24);
    plt.ylabel(r'$\eta$' + '    ' + r'$Pa \cdot s$', fontsize=24);
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
    
    print(parameters)
    
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
    plt.xlabel(r'$\gamma$' + '    ' + r'$1/s$', fontsize=24);
    plt.ylabel(r'$N1$' + '    ' + r'$Pa$', fontsize=24);
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
    if len(y_str) > 1:
        x_strlabel = '_'.join(x_str[0].split('_')[0:-1])
        y_strlabel = '_'.join(y_str[0].split('_')[0:-1])
    else:
        x_strlabel = x_str
        y_strlabel = y_str
    ax0.set_xlabel(str(x_strlabel) + '$    ' + str(x_units), fontsize=24);
    ax0.set_ylabel(str(y_strlabel) + '$    ' + str(y_units), fontsize=24);

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