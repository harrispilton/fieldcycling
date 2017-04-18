#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from bokeh.charts import Scatter, output_file, show
from bokeh.sampledata.autompg import autompg as df
from bokeh.models.widgets import Panel, Tabs
from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Dropdown
from tornado.ioloop import IOLoop
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, Slider
from bokeh.server.server import Server
import stelardatafile as sdf

path=os.path.join(os.path.curdir,'data')
polymer=sdf.StelarDataFile('297K.sdf',path)
polymer.sdfimport()
nr_experiments = polymer.get_number_of_experiments()

io_loop = IOLoop.current()
ie = 1

def modify_doc(doc):
    parameters=polymer.getparameter(ie)
    bs=int(parameters['BS'])
    try:
        nblk=int(parameters['NBLK'])
    except:
        nblk=1;
    ns=int(parameters['NS'])
    try:
        dw=parameters['DW']*1e-6 #dwell time is in [mu s]
    except:
        dw=1
    temperature=parameters['TEMP']
    fid=pd.DataFrame(polymer.getdata(ie),index=np.linspace(dw,dw*bs*nblk,bs*                 nblk), columns=['real', 'im'])/ns
    fid['magnitude']=( fid['real']**2 + fid['im']**2 )**0.5 # last two lines may represent a seperate method
    tau=np.logspace(-3,np.log10(5*parameters['T1MX']),nblk) #as a dummy
    startpoint=int(0.05*bs)-1
    endpoint=int(0.1*bs)
    phi=np.zeros(nblk)
    for blk in range(nblk):
        start=startpoint + blk * bs
        end=endpoint + blk * bs
        phi[blk]=fid['magnitude'].iloc[start:end].sum() / (endpoint-startpoint)
    df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])

    source_fid = ColumnDataSource(data=ColumnDataSource.from_df(fid))
    source_df = ColumnDataSource(data=ColumnDataSource.from_df(df))

    
    p1 = figure(plot_width=300, plot_height=300)
    p1.line('index', 'im', source=source_fid, color='blue')
    p1.line('index', 'real', source=source_fid, color='green')
    p1.line('index', 'magnitude', source=source_fid, color='red')


    p2 = figure(plot_width=300, plot_height=300)
    p2.circle_cross('tau', 'phi', source=source_df, color="navy")
        
    menu = [("Experiment nr {:4d}".format(ie), "{:4d}".format(ie)) for ie in range(1, nr_experiments+1)]
    dropdown = Dropdown(label="choose experiment number", menu=menu)    

    def callback(attr, old, new):
        ie = int(new)
        parameters=polymer.getparameter(ie)
        bs=int(parameters['BS'])
        try:
            nblk=int(parameters['NBLK'])
        except:
            nblk=1;
        ns=int(parameters['NS'])
        try:
            dw=parameters['DW']*1e-6 #dwell time is in [mu s]
        except:
            dw=1
        temperature=parameters['TEMP']
        fid=pd.DataFrame(polymer.getdata(ie),index=np.linspace(dw,dw*bs*nblk,bs*                 nblk), columns=['real', 'im'])/ns
        fid['magnitude']=( fid['real']**2 + fid['im']**2 )**0.5 # last two lines may represent a seperate method
        tau=np.logspace(-3,np.log10(5*parameters['T1MX']),nblk) #as a dummy
        startpoint=int(0.05*bs)-1
        endpoint=int(0.1*bs)
        phi=np.zeros(nblk)  
        for blk in range(nblk):
            start=startpoint + blk * bs
            end=endpoint + blk * bs
            phi[blk]=fid['magnitude'].iloc[start:end].sum() / (endpoint-startpoint)
        df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])

        source_fid.data = ColumnDataSource.from_df(fid)
        source_df.data = ColumnDataSource.from_df(df)

    dropdown.on_change('value', callback)
    doc.add_root(column(dropdown, p1, p2))

bokeh_app = Application(FunctionHandler(modify_doc))

server = Server({'/': bokeh_app}, io_loop=io_loop)
server.start()

if __name__ == '__main__':
    print('Opening Bokeh application on http://localhost:5006/')

    io_loop.add_callback(server.show, "/")
io_loop.start()