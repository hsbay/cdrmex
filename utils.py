# CDRMEX/utilspy
# utility functions for CDRMEx, CC-BY-40 2021, @safiume

# estimate CS routine heavily borrows from pymagicccorepy
# pymagicc, AGPL-30, 2021

import sys, re
import shutil
import subprocess
import warnings
from copy import deepcopy
from os import listdir, makedirs
from os.path import abspath, basename, dirname, exists, isfile, join
from subprocess import PIPE
from tempfile import mkdtemp
import numpy as np
import expectexception
from datetime import datetime
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
# %matplotlib inline
import f90nml
from dateutil.relativedelta import relativedelta
from openscm_units import unit_registry
from scmdata import run_append

import pymagicc
from pymagicc.config import _wine_installed, config
from pymagicc.io import MAGICCData, read_cfg_file
from pymagicc.io.utils import _get_openscm_var_from_filepath
from pymagicc.scenarios import zero_emissions
from pymagicc.utils import get_date_time_string

# Plothelper, Set up matplotlib defs
# plthelpr(Plot axes, plot, setables='foo')

def gline(profile):
    # chose: full, pub
    # Graph start and end, adjust the dates below
    if profile == 'full':
        start = 1700
        end = 2700
    elif profile == 'pub':
        start = 1850
        end = 2125
    elif profile == 'emiss':
        start = 2010
        end = 2125
    else:
        start = graphstart
        end = graphend
    x = datetime(start,1,1,0), datetime(end,1,1,1)
    return(x)

def plthelpr(pltax,plt,**kwargs):
    x = gline(kwargs['profile'])
    mlocator = mdates.YearLocator(50, month=1, day=1)
    minloc = mdates.YearLocator(10, month=1, day=1)
    formatter = mdates.ConciseDateFormatter(mlocator)
    pltax.set_title(var)                        
    pltax.set_xlim(x)
    pltax.set_ylabel('('+ylab+')')
    pltax.tick_params(which='major', width=0.75, length=6)
    pltax.xaxis.set_major_locator(mlocator)
    pltax.xaxis.set_major_formatter(formatter)
    pltax.xaxis.set_minor_locator(minloc)
    if kwargs['profile'] == 'pub' or kwargs['profile'] == 'emiss':
        wh = 'both'
    else:
        wh = 'major'

    if 'clr' in kwargs:
        clr = kwargs['clr']
    else:
        clr = '.80'

    pltax.grid(which=wh, linewidth=0.5, color=clr)
    plt.tight_layout(pad=0.6, h_pad=1.5)

tstcfg = {
        'co2_switchfromconc2emis_year' : 30000,
        'rf_tropoz_constantafteryr' : 5000,
        'rf_stratoz_constantafteryr' : 5000,
        'rf_mhalo_constantafteryr' : 5000,
        'rf_fgas_constantafteryr' : 5000,
        'rf_landuse_constantafteryr' : 5000,
        'rf_mineraldust_constantafteryr' : 5000,
        'n2o_switchfromconc2emis_year': 5000,
        'fgas_switchfromconc2emis_year': 5000,
        'mhalo_switch_conc2emis_yr' : 5000,
        'rf_total_runmodus':'CO2',
        'core_heatxchange_landocean':1,
        'out_emissions': 1,
        'out_parameters' : 1,
        'out_inverseemis' : 1,
}

def diagnose_tcr_ecs_tcre(direction, **kwargs):

    # more generic handling of positive and negative ECS testing
    # borrows heavily from pymagicc and carries AGPL3 license.

    # diagnose_tcr_ecs_tcre([pos|neg], **kwargs):

    ecscfg = { 'startyear' : 1795,
        'endyear' : 4321,
        'core_climatesensitivity' : 3.6, }
    ecscfg['core_climatesensitivity'] = kwargs['CORE_CLIMATESENSITIVITY']
    ecscfg['core_delq2xco2'] = kwargs['CORE_DELQ2XCO2']
    ecscfg['rf_total_constantafteryr'] = 2000

    tcrcfg = { 'startyear' : 1750,
        'endyear' : 2573,
        'core_climatesensitivity' : 3.6, }
    tcrcfg['core_climatesensitivity'] = kwargs['CORE_CLIMATESENSITIVITY']
    tcrcfg['rf_total_constantafteryr'] = 3000

    ecs_res = diagnose_ecs(direction, **ecscfg,**tstcfg)
    tcr_tcre_res = diagnose_tcr_tcre(direction, **tcrcfg,**tstcfg)

    out = {**ecs_res, **tcr_tcre_res}
    out['timeseries'] = run_append(
      [ecs_res['timeseries'], tcr_tcre_res['timeseries'],]
    )

    return out

def diagnose_ecs(direction, **kwargs):
    posecstest = { 'file_co2_conc' : 'ABRUPT2XCO2_CO2_CONC.IN',
       'testscen': 'abrupt-2xCO2' }
    negecstest = { 'file_co2_conc' : 'ABRUPT0P5XCO2_CO2_CONC.IN',
       'testscen' : 'abrupt-0p5xCO2' }

    if 'pos' in direction:
       selectedtst = posecstest
    elif 'neg' in direction:
       selectedtst = negecstest

    print('Calculating ECS from {}.'.format(selectedtst['testscen']))

    timeseries = pymagicc.run(
        scenario=None,
        only=[
            'Atmospheric Concentrations|CO2',
            'INVERSEEMIS',
            'Radiative Forcing',
            'Surface Temperature',
        ], **kwargs, file_co2_conc = selectedtst['file_co2_conc']
    )
    # drop all the irrelevant inverse emissions
    timeseries = timeseries.filter(
        variable='Inverse Emissions*', level=1, keep=False
    )
    timeseries['scenario'] = selectedtst['testscen']

    ecs = get_ecs_from_diagnosis_results(timeseries)
    return {'ecs': ecs, 'timeseries': timeseries}

def diagnose_tcr_tcre(direction, **kwargs):
    postcrtest = { 'file_co2_conc' : '1PCTCO2_CO2_CONC.IN',
       'testscen' : '1pctCO2' }
    negtcrtest = { 'file_co2_conc' : '1PCTCDR_CO2_CONC.IN',
       'testscen' : '1pctCO2-cdr' }

    if 'pos' in direction:
       selectedtst = postcrtest
    elif 'neg' in direction:
       selectedtst = negtcrtest

    print('Calculating TCR & TCRE from {}.'.format(selectedtst['testscen']))

    timeseries = pymagicc.run(
        scenario=None,
        only=[
            'Atmospheric Concentrations|CO2',
            'INVERSEEMIS',
            'Radiative Forcing',
            'Surface Temperature',
        ], **kwargs, file_co2_conc = selectedtst['file_co2_conc']
    )
    # drop all the irrelevant inverse emissions
    timeseries = timeseries.filter(
        variable='Inverse Emissions*', level=1, keep=False
    )
    timeseries['scenario'] = selectedtst['testscen']
    tcr, tcre = get_tcr_tcre_from_diagnosis_results(timeseries)

    return {'tcr': tcr, 'tcre': tcre, 'timeseries': timeseries}

def get_ecs_from_diagnosis_results( results_ecs_run):
    global_co2_concs = results_ecs_run.filter(
        variable='Atmospheric Concentrations|CO2', region='World'
    )
    ecs_time, ecs_start_time = get_ecs_ecs_start_yr_from_CO2_concs(
        global_co2_concs
    )

    global_total_rf = results_ecs_run.filter(
        variable='Radiative Forcing', region='World'
    )

    global_temp = results_ecs_run.filter(
        variable='Surface Temperature', region='World'
    )

    ecs = float(global_temp.filter(time=ecs_time).values.squeeze())
    ecs = abs(ecs)
    unit = global_temp.get_unique_meta('unit', no_duplicates=True)
    ecs = ecs * unit_registry(unit)

    return ecs

def get_tcr_tcre_from_diagnosis_results( results_tcr_tcre_run):
    global_co2_concs = results_tcr_tcre_run.filter(
            variable='Atmospheric Concentrations|CO2', region='World'
        )
    (tcr_time, tcr_start_time) = get_tcr_tcr_start_yr_from_CO2_concs(
            global_co2_concs
        )

    global_inverse_co2_emms = results_tcr_tcre_run.filter(
            variable='Inverse Emissions|CO2|MAGICC Fossil and Industrial',
            region='World',
    )

    global_total_rf = results_tcr_tcre_run.filter(
            variable='Radiative Forcing', region='World'
    )

    global_temp = results_tcr_tcre_run.filter(
            variable='Surface Temperature', region='World'
    )

    tcr = float(global_temp.filter(time=tcr_time).values.squeeze())
    tcr_unit = global_temp.get_unique_meta('unit', no_duplicates=True)
    tcr = tcr * unit_registry(tcr_unit)

    tcre_cumulative_emms = float(
        global_inverse_co2_emms.filter(
         year=range(tcr_start_time.year, tcr_time.year)
        ).values.sum()
    )
    emms_unit = global_inverse_co2_emms.get_unique_meta('unit', no_duplicates=True)
    years = global_inverse_co2_emms['year'].values.squeeze()
    tcre_cumulative_emms_unit = unit_registry(emms_unit) * unit_registry('yr')
    tcre_cumulative_emms = tcre_cumulative_emms * tcre_cumulative_emms_unit
    
    tcre = tcr / tcre_cumulative_emms

    return tcr, tcre

def get_ecs_ecs_start_yr_from_CO2_concs( df_co2_concs):
    co2_concs = df_co2_concs.timeseries()
    co2_conc_0 = co2_concs.iloc[0, 0]
    t_start = co2_concs.columns.min()
    t_end = co2_concs.columns.max()

    ecs_start_time = co2_concs.iloc[
              :, co2_concs.values.squeeze() == co2_conc_0
    ].columns[-1]

    ecs_time = df_co2_concs['time'].iloc[-1]

    return ecs_time, ecs_start_time

def get_tcr_tcr_start_yr_from_CO2_concs( df_co2_concs):
    co2_concs = df_co2_concs.timeseries()
    co2_conc_0 = co2_concs.iloc[0, 0]
    t_start = co2_concs.columns.min()
    t_end = co2_concs.columns.max()

    tcr_start_time = co2_concs.iloc[
        :, co2_concs.values.squeeze() > co2_conc_0
    ].columns[0] - relativedelta(years=1)
    tcr_time = tcr_start_time + relativedelta(years=70)

    return tcr_time, tcr_start_time


