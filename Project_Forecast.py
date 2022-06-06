#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 15:34:20 2022

@author: maximilienpeters
"""
import numpy as np
import matplotlib.pyplot as plt
import datetime
from dateutil.relativedelta import relativedelta
import matplotlib.dates as mdates
from scipy.stats import norm
from matplotlib import rc
import matplotlib.animation as animation

plt.style.use('seaborn-whitegrid')

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica'], 'size': 12})
rc('text', usetex=True)

#------------------------------------------------------------------------------------------

file = open('time_series.txt', 'r') #see data file "Time_series.txt" for data format

list_of_lists = []
for line in file:
    stripped_line = line.strip()
    line_list = stripped_line.split()
    list_of_lists.append(line_list)

#------------------------------------------------------------------------------------------

def Time_series(start_year, end_year, m):#get a time series of the desired month from the initial time series

    dates_list = []
    Delta_years = (end_year+1-start_year)
    
    if end_year == 2022 and m > 5:
        Delta_years = Delta_years-1
    
    SIE = np.zeros(Delta_years)
    dates = np.zeros(Delta_years)
    years_count = 0

    for i in range(12*Delta_years):
            
        if (i-(m-1)) % 12 == 0 and (i-(m-1))>=0: 
            SIE[years_count] = int(list_of_lists[i][3])
            dates_list = dates_list + [datetime.datetime(int(list_of_lists[i][0]),int(list_of_lists[i][1]),int(list_of_lists[i][2]))]
            years_count = years_count + 1
                            
    SIE = np.where(SIE == -999, np.nan, SIE)
    dates = np.array(dates_list)
    
    return(dates, SIE)

def Time_series_total():
    
    dates_list = []
    
    SIE = np.zeros(len(list_of_lists))
    dates = np.zeros(len(list_of_lists))

    for i in range(len(list_of_lists)):
            
        SIE[i] = int(list_of_lists[i][3])
        dates_list = dates_list + [datetime.datetime(int(list_of_lists[i][0]),int(list_of_lists[i][1]),int(list_of_lists[i][2]))]
                            
    SIE = np.where(SIE == -999, np.nan, SIE)
    dates = np.array(dates_list)
    
    return(dates, SIE)


def plot_time_series_total():
    
    fig, ax = plt.subplots(constrained_layout=True)
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)  
    
    ax.plot(dates, SIE)
    
    ax.set_xlabel('Year')
    ax.set_ylabel('Sea ice extent (km$^2$)')
    
    plt.savefig("SIE_time_series_total", dpi = 900)
    plt.show()


def trend(start_year, end_year, SIE):
    
    a_num = 0
    a_den = 0
    Delta_years = (end_year+1-start_year)
    x_mean = start_year + Delta_years/2 - 0.5
    x = np.linspace(start_year, end_year, Delta_years)
    SIE_sum = 0
    n = start_year - 1979
    
    for j in range(Delta_years-n):
        SIE_sum = SIE_sum + SIE[j+n]
        
    SIE_mean = SIE_sum/(Delta_years)
    
    for i in range(Delta_years-n):
        a_num = a_num + ((x[i]-x_mean)*(SIE[i+n]-SIE_mean))
        a_den = a_den + ((x[i]-x_mean)*(x[i]-x_mean))
    
    a = a_num/a_den
    b = SIE_mean - a*x_mean
    Y = b + a*x
    Trend_last_value = Y[-1] + a 
    
    return(Y, Trend_last_value)
    

def Mean_Var(start_year, end_year, month):
    
        
    Delta_years = (end_year+1-start_year)
    
    if month == 9:
        dates_month, SIE_month = Time_series(start_year, end_year-1, month)
    else:
        dates_month, SIE_month = Time_series(start_year, end_year, month)
        
    SIE_month_sums = np.zeros(Delta_years-2)
    SIE_month_means = np.zeros(Delta_years-2)
    VAR_SIE_month = np.zeros(Delta_years-2)
    
    if month > 9:
        SIE_month = np.insert(SIE_month[:Delta_years-1], 0, np.nan) #Add NaN for October-November-December of 1978 explained by the fact that we don't have data for 1978
    
    for i in range(Delta_years-2):
        n = 0
        m = 0
        
        for j in range(i+2):

            if not np.isnan(SIE_month[j]):
                SIE_month_sums[i] = SIE_month_sums[i] + SIE_month[j]
                n = n + 1
            
        SIE_month_means[i] = SIE_month_sums[i]/n
        
        for k in range(i+2):
            if not np.isnan(SIE_month[k]):
                VAR_SIE_month[i] = VAR_SIE_month[i] + (SIE_month[k] - SIE_month_means[i])**2
                m = m + 1
        
        VAR_SIE_month[i] = VAR_SIE_month[i]/m
    
    return(SIE_month_means, VAR_SIE_month, SIE_month)


def Forecast_system(start_year, end_year, n_month): 
    
    Delta_years = (end_year+1-start_year)
    Standard_dev = np.zeros(Delta_years-2)
    Mean_values = np.zeros(Delta_years-2)
    Anomaly = np.zeros(Delta_years-2)
    
    SIE_month_means = np.zeros((n_month,Delta_years-2))
    VAR_SIE_month = np.zeros((n_month,Delta_years-2))
    SIE_month = np.zeros((n_month,Delta_years))
    
    for j in range(n_month):
        
        month = 5 - (n_month-1) + j
        if month <= 0:
            month = 12 + month
        
        SIE_month_means[j], VAR_SIE_month[j], SIE_month[j] = Mean_Var(start_year, end_year, month)
    
    SIE_sept_means, VAR_SIE_sept, SIE_sept = Mean_Var(start_year, end_year, 9)
    
    for i in range(Delta_years-2):
    
        n = 0
        if n_month > 1:
            Standard_dev[i] = ((VAR_SIE_sept[i] + np.sum(VAR_SIE_month.T[i]))/(i+2))**(1/2)#The .T transpose the matrix of month in order to sum correclty the variance of each month until the year i
        else: 
            Standard_dev[i] = ((VAR_SIE_sept[i] + VAR_SIE_month[0][i])/(i+2))**(1/2)
        
        for j in range(n_month):
            if not np.isnan(SIE_month[j][i+2]):
                Anomaly[i] = Anomaly[i] + (SIE_month[j][i+2] - SIE_month_means[j][i])
                n = n + 1
        
        Anomaly[i] = Anomaly[i]/n #warning for n_month = 1 or 2 is expected -> artifical value of mean_value is automaticly assigned for i = 5
        Mean_values[i] = SIE_sept_means[i] + Anomaly[i]

    
    if n_month == 1 or n_month == 2:
        Mean_values[5] = (Mean_values[4] + Mean_values[6])/2
        
    return(Mean_values, Standard_dev, VAR_SIE_sept)


def Forecast_correction(data_start_year, data_end_year):
    
    if data_end_year == 2022:
        bias = mu_1[:-1] - SIE_sept[2:]
    else:
        bias = mu_1 - SIE_sept[2:]
        
    mean_bias = np.sum(bias)/np.size(bias) #first correction bias of the mean
    
    if data_end_year == 2022:
        Trend_obs, Trend_last_value = trend(data_start_year, data_end_year-1, SIE_sept)
        Trend_obs = np.append(Trend_obs, Trend_last_value)
    else:
        Trend_obs, Trend_last_value = trend(data_start_year, data_end_year, SIE_sept)
    
    new_mu_1 = (mu_1 - mean_bias) #bias of the mean
    
    Trend_mean_values, E_1 = trend(data_start_year+2, data_end_year, new_mu_1)
    new_mu_3 = mu_1 - Trend_mean_values + Trend_obs[2:] #bias of the trend
    
    if data_end_year == 2022:
        bias = new_mu_3[:-1] - SIE_sept[2:]
    else:
        bias = new_mu_3 - SIE_sept[2:]
    
    mean_bias = np.sum(bias)/np.size(bias)
    new_mu_4 = (new_mu_3 - mean_bias) #last correction bias of the mean
    
    return(new_mu_4)


def plot_Bias():
    
    fig4, ax4 = plt.subplots(constrained_layout=True)
    
    x = np.linspace(3e6, 9e6, 100)
    
    ax4.plot(mu, SIE_sept[2:], ".")
    ax4.plot(x, x, "--")
    ax4.set_title('Bias')
    plt.savefig("Bias", dpi = 900)
    plt.show()


def Event_forecast(start_year, end_year, Event):#end_year = last verification year / Event = 1 -> trend line event, Event = 3 -> less than previous year event
    
    Delta_years = (end_year+1-start_year)
    E = np.zeros(Delta_years-2)
    Event_occ = np.zeros(Delta_years-2)
    Oberved_frequencies = np.zeros(Delta_years-2)
    P_event = np.zeros(Delta_years-2)
    
    if Event == 1:
    
        for i in range(Delta_years-2):
        
            Y,E[i] = trend(start_year, start_year+1+i, SIE_sept)
        
            if SIE_sept[i+2] < E[i]: 
                Event_occ[i] = 1
            else:
                Event_occ[i] = 0
        
            P_event[i] = norm.cdf((E[i] - mu[i])/sigma[i])
        
        for j in range(Delta_years-2):
        
            for k in range(j+1):
            
                Oberved_frequencies[j] = Oberved_frequencies[j] + Event_occ[k]
        
            Oberved_frequencies[j] = Oberved_frequencies[j]/(j+1)
    
    elif Event == 3:
        
        for i in range(Delta_years-2):
            
            if SIE_sept[i+2] < SIE_sept[i+1]: 
                Event_occ[i] = 1
            else:
                Event_occ[i] = 0
            
            P_event[i] = norm.cdf((SIE_sept[i+1] - mu[i])/sigma[i])
            
        for j in range(Delta_years-2):
        
            for k in range(j+1):
            
                Oberved_frequencies[j] = Oberved_frequencies[j] + Event_occ[k]
        
            Oberved_frequencies[j] = Oberved_frequencies[j]/(j+1)
        
    return(P_event, Event_occ, Oberved_frequencies)


def Proba_forcast(start_year, end_year):
    
    fig3, ax3 = plt.subplots()
    
    years = np.arange(start_year,end_year+1)
    
    P_prediction = P
    
    if end_year == 2022:
        P_prediction = np.append(P, norm.cdf((SIE_sept[-1] - mu[-1])/sigma[-1]))
    
    P_predi_E = np.zeros(np.size(P_prediction))
    P_predi_N = np.zeros(np.size(P_prediction))
    P_predi_L = np.zeros(np.size(P_prediction))
    
    for i in range(np.size(P_prediction)):
        if i == np.size(P_prediction)-1 and end_year == 2022:
            P_predi_L[i] = P_prediction[i]
        elif E_O[i] == 1:
            P_predi_E[i] = P_prediction[i]
        else:
            P_predi_N[i] = P_prediction[i]
            
    ax3.bar(years, P_predi_E, color = "#50D51F", label = 'Event')
    ax3.bar(years, P_predi_N, color = "#FF5733", label = 'No event')
    ax3.bar(years, P_predi_L, color = "#1F92D5", label = '2022 Prediction')
        
    ax3.set_xlabel('Year')
    ax3.set_ylabel('Predicted probability')
    plt.legend(frameon = True, facecolor = 'white', framealpha = 0.7, loc = "upper left")
    plt.savefig("Predi_proba", dpi = 900)
    plt.show()
    
    return(P_prediction[-1])


def Brier_skill_score():
    
    BS = 0
    BS_ref = 0
    
    for i in range(np.size(P)):
        BS = BS + (P[i] - E_O[i])**2
        
    BS_ref = (O_f[-1]) * (1 - (O_f[-1]))
    BS = BS/np.size(P)
    
    BSS = 1 - BS/BS_ref

    return(BSS, BS, BS_ref)


def plot_reliability_diagram():
    
    fig2, ax2 = plt.subplots(constrained_layout=True)
    
    ax2.plot(P, O_f, ".")
    ax2.set_title('Reliability diagram 1979-2021')
    plt.savefig("Reliability diagram", dpi = 900)
    plt.show()
    

def plot_time_series(trend_start_year, trend_end_year):
    
    fig, ax = plt.subplots(constrained_layout=True)
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    
    ax.plot(dates_sept, SIE_sept, label = 'Observed SIE')
    ax.plot(dates_mu, mu, label = 'Predicted SIE')
    ax.fill_between(dates_mu, mu - 2*sigma, mu + 2*sigma, alpha = 0.5, label = "95\% CI")
    
    Y, E = trend(trend_start_year, trend_end_year, SIE_sept)
    
    ax.plot(dates_sept[trend_start_year-1979:trend_end_year+1-1979], Y, label = 'Trend of Observations')
    ax.set_xlabel('Year')
    ax.set_ylabel('Sea ice extent (km$^2$)')
    plt.legend()
    plt.savefig("SIE_time_series", dpi = 900)
    plt.show()

#------------------------------------------------------------------------------------------

start_year = 1979
end_year = 2022

Total_number_of_month = 8
Event_number = 3

dates, SIE = Time_series_total()

plot_time_series_total()

if end_year <= 2021:
    dates_sept, SIE_sept = Time_series(start_year, end_year, 9)
    dates_mu = dates_sept[2:]
elif end_year == 2022:
    dates_sept, SIE_sept = Time_series(start_year, end_year-1, 9)
    dates_mu = np.append(dates_sept[2:], dates_sept[-1] + relativedelta(years=1))

mu_1, sigma, VAR_sept = Forecast_system(start_year, end_year, Total_number_of_month)

mu = Forecast_correction(start_year, end_year)


if end_year <= 2021:
    plot_time_series(start_year,end_year)
elif end_year == 2022:
    plot_time_series(start_year,end_year-1)

#plot_Bias()

if end_year <= 2021:
    P, E_O, O_f = Event_forecast(start_year, end_year, Event_number) 
elif end_year == 2022:
    P, E_O, O_f = Event_forecast(start_year, end_year-1, Event_number) 

#plot_reliability_diagram()

P_event_2022 = Proba_forcast(start_year+2, end_year)

BSS, BS, BS_ref = Brier_skill_score()

print("Predicted Probability of event in 2022 is ", P_event_2022)

print("Predicted September SIE for 2022 is ", mu[-1], " km2")

print("BSS = ", BSS, "\n", "BS = ", BS, "\n", "BS_ref = ", BS_ref)

#------------------------------------------------------------------------------------------
def Video_time_series(): 
    array_1 = np.append(SIE_sept, np.nan)
    array_2 = np.insert(np.insert(mu, 0, np.nan), 0, np.nan)
    sigma2 = np.insert(np.insert(sigma, 0, np.nan), 0, np.nan)
    array_3 = array_2 - 2*sigma2
    array_4 = array_2 + 2*sigma2

    dates = np.append(dates_sept, dates_sept[-1] + relativedelta(years=1))
    
    array_y1, array_y2, array_y3, array_y4, array_x = [], [], [], [], []
    

    def update(n):
        array_y1.append(array_1[n])
        array_y2.append(array_2[n])
        array_y3.append(array_3[n])
        array_y4.append(array_4[n])
        array_x.append(dates[n])
        ax.fill_between(array_x, array_y3, array_y4, color = "#98D0F0")
        ax.plot(array_x, array_y1, color = "tab:blue", label = 'Observed SIE')
        ax.plot(array_x, array_y2, color = "tab:orange", label = 'Predicted SIE')
        ax.set_xlabel('Year')
        ax.set_ylabel('Sea ice extent (km$^2$)')

    fig, ax = plt.subplots(constrained_layout=True)
    locator = mdates.AutoDateLocator()
    formatter = mdates.ConciseDateFormatter(locator)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_ylim(4e6, 9e6)
    ax.set_xlim(dates[0] - relativedelta(years=2), dates[-1] + relativedelta(years=2))
    
    ani = animation.FuncAnimation(fig, update, frames = 44, interval=400)
    writervideo = animation.FFMpegWriter(fps=8)
    ani.save('SIE_predi_movie.mp4', writer=writervideo, dpi = 400)
    plt.close()

Video_time_series()


