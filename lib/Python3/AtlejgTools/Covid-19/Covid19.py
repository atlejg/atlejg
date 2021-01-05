#!/usr/bin/env python
# coding: utf-8

# # 1. Initialise

# In[1]:


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
# get_ipython().run_line_magic('matplotlib', 'inline')

from covid.api import CovId19Data


# In[304]:


# pip install covid-data-api
# https://pypi.org/project/covid-data-api/

api = CovId19Data(force=True)


# # 2. Read country data sets

# In[305]:


def get_history(Country):
    res = api.get_history_by_country(Country)
    data = res[Country]
    data = data['history']
    covid = pd.DataFrame.from_dict(data, orient='index')
    covid.index = pd.to_datetime(covid.index)
    covid.replace(r'na', np.nan,inplace=True)
    covid.change_confirmed = pd.to_numeric(covid.change_confirmed, downcast='float')
    covid.change_deaths = pd.to_numeric(covid.change_deaths, downcast='float')
    covid['time_to_double'] = np.log(2)/np.log(1 + covid.change_confirmed)
    return covid


# In[306]:


def get_south_korea():
    res = api.get_history_by_country('korea, south')
    data = res['korea_south']
    data = data['history']
    covid = pd.DataFrame.from_dict(data, orient='index')
    covid.index = pd.to_datetime(covid.index)
    covid.replace(r'na', np.nan,inplace=True)
    covid.change_confirmed = pd.to_numeric(covid.change_confirmed, downcast='float')
    covid.change_deaths = pd.to_numeric(covid.change_deaths, downcast='float')
    covid['time_to_double'] = np.log(2)/np.log(1 + covid.change_confirmed)
    return covid


# In[307]:


covid_norge = get_history('norway')
covid_sveits = get_history('switzerland')
covid_italy = get_history('italy')
covid_usa = get_history('us')
covid_singapore = get_history('singapore')
covid_southkorea = get_south_korea()


# # 3. Plots of Growth Rates

# In[314]:


kwargs = dict(xlim=('2020-03-02','2020-03-30'),ylim=(0,1))
ax = covid_norge.change_confirmed.plot(figsize=(14,5), lw=3, **kwargs)
ax.set_ylabel('Daily change of Confirmed Cases')
covid_sveits.change_confirmed.plot(ax=ax, **kwargs)
covid_italy.change_confirmed.plot(ax=ax, **kwargs)
covid_usa.change_confirmed.plot(ax=ax, **kwargs)
# covid_sverige.change_confirmed.plot(ax=ax, **kwargs)
ax.legend(['Norway','Switzerland','Italy','USA'])
plt.show()


# In[313]:


kwargs = dict(xlim=('2020-03-02','2020-03-30'),ylim=(12,0))
ax = covid_norge.time_to_double.rolling(3).mean().plot(figsize=(14,5), lw=3, **kwargs)
ax.set_ylabel('Doubling time Confirmed Cases, days')
covid_sveits.time_to_double.rolling(3).mean().plot(ax=ax, **kwargs)
covid_italy.time_to_double.rolling(3).mean().plot(ax=ax, **kwargs)
covid_usa.time_to_double.rolling(3).mean().plot(ax=ax, **kwargs)
ax.legend(['Norway','Switzerland','Italy','USA']);
plt.title('Doubling time (rolling 3-day average)')
ax.text(0.05, 0.2, 'less than 4 days: fast \n4-8 days: moderate, not slowing down \n8-30 days: slowing down', transform=ax.transAxes)
ax.get_figure().patch.set_facecolor('w')
plt.show()


# Daily change of confirmed cases:   $c$  
#   
# Doubling time  $x$ in days:  $2 = (1 + c)^x$  $$x = \frac{ln(2)}{ln(1 + c)}$$
# 
# **Doubling time of cases in days:**  
# less than 4 days: fast  
# 4-8 days: moderate, not slowing down  
# 8-30 days: slowing down  
# \>30 days: slowing down a lot  

# In[312]:


kwargs = dict(ylim=(12,0))
ax = covid_norge['2020-03-06':].reset_index().time_to_double.rolling(3).mean().plot(figsize=(14,5), lw=3, **kwargs)
covid_sveits['2020-03-05':].reset_index().time_to_double.rolling(3).mean().plot(ax=ax, **kwargs)
covid_italy['2020-02-23':].reset_index().time_to_double.rolling(3).mean().plot(ax=ax, **kwargs)
covid_usa['2020-03-02':].reset_index().time_to_double.rolling(3).mean().plot(ax=ax, **kwargs)
#covid_singapore['2020-02-29':].reset_index().time_to_double.plot(ax=ax, **kwargs)
ax.set_ylabel('Doubling time Confirmed Cases, days')
ax.set_xlabel('Days since the total confirmed cases of COVID-19 reached 100')
ax.legend(['Norway','Switzerland','Italy','USA'])
plt.show()


# In[141]:


def f(days, double):
    return 100*2**(days/double)


# In[200]:


double10 = pd.DataFrame({'days':[0, 10, 30, 50]})
double10['cases'] = f(double10.days,10)
double10.set_index('days',inplace=True)
double5 = pd.DataFrame({'days':[0, 10, 25, 50]})
double5['cases'] = f(double5.days,5)
double5.set_index('days',inplace=True)
double2 = pd.DataFrame({'days':[0, 10, 50]})
double2['cases'] = f(double2.days,2)
double2.set_index('days',inplace=True)


# In[311]:


kwargs = dict(logy=True,ylim=(100,100000),xlim=(0,40))
ax = covid_norge['2020-03-06':].reset_index().confirmed.plot(figsize=(10,7), lw=3, **kwargs)
covid_sveits['2020-03-05':].reset_index().confirmed.plot(ax=ax, **kwargs)
covid_italy['2020-02-22':].reset_index().confirmed.plot(ax=ax, **kwargs)
covid_usa['2020-03-02':].reset_index().confirmed.plot(ax=ax, **kwargs)
covid_singapore['2020-02-29':].reset_index().confirmed.plot(ax=ax, **kwargs)
covid_southkorea['2020-02-20':].reset_index().confirmed.plot(ax=ax, **kwargs)
double10.cases.plot(ax=ax, style='k--', alpha=0.2, **kwargs)
double5.cases.plot(ax=ax, style='k--', alpha=0.2, **kwargs)
double2.cases.plot(ax=ax, style='k--', alpha=0.2, **kwargs)
ax.set_ylabel('Total confirmed cases of COVID-19')
ax.set_xlabel('Days since the total confirmed cases of COVID-19 reached 100')
ax.legend(['Norway','Switzerland','Italy','USA','Singapore','South Korea']);
plt.text(30, 1000, 'doubling every 10 days', rotation=15, alpha=0.2, rotation_mode='anchor')
plt.text(25, 4000, 'doubling every 5 days', rotation=30, alpha=0.2, rotation_mode='anchor')
plt.text(10, 4000, 'doubling every 2 days', rotation=55, alpha=0.2, rotation_mode='anchor')
ax.annotate('March 12$^{th}$', xy=(6,702), xycoords='data',xytext=(3,4000),textcoords='data', arrowprops=dict(arrowstyle='->',color='b',alpha=0.4),color='b',alpha=0.4)
ax.get_figure().patch.set_facecolor('w')
plt.show()


# In[ ]:


# fig.patch.set_facecolor('w')


# In[ ]:





# In[292]:


# covid_southkorea['2020-02-01':]


# In[278]:


api.show_available_countries()


# In[288]:


api.get_history_by_country('korea, south')


# In[ ]:




