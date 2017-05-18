#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os
import os.path as osp
import pandas as pd
import numpy as np
from bokeh.models.widgets import Panel, Tabs
from bokeh.plotting import figure, ColumnDataSource
from bokeh.palettes import Spectral5, viridis  # @UnresolvedImport
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Select, Button, DataTable, TableColumn, Slider, RangeSlider, TextInput#, ColumnDataSource
from bokeh.models.callbacks import CustomJS
from tornado.ioloop import IOLoop
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
from bokeh.layouts import column, row
#from bokeh.models import ColumnDataSource#, Slider, RangeSlider, TextInput
from bokeh.server.server import Server
from scipy.optimize import leastsq
import logging

from utils import model_exp_dec, fun_exp_dec, get_mag_amplitude, magnetization_fit
import stelardatafile as sdf

# Set up a logger.
logger = logging.getLogger(__name__)
#set this to loggin.INFO if the console gets to crowded
logger.setLevel(logging.DEBUG)
# create a file handler
handler = logging.FileHandler('logfile.log')
handler.setLevel(logging.INFO)
# create a logging format
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(funcName)s - %(message)s')
handler.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(handler)

# in the plot 4 use following impo
SIZES = list(range(6, 22, 3)) # for some sizes
COLORS = Spectral5 # for some colors (for more colors use viridis(integer))

def load_data(sdf_file, path=None):
    # TODO: give path as a parameter to use this function for other datafiles
    #specify and import data file
    if not path:
        path=os.path.join(os.path.curdir,'data') # default path, if it's not set as an argument
    # initialize StelarDataFile Object
    polymer=sdf.StelarDataFile(sdf_file, path, logger=logger)
    polymer.sdfimport()
    nr_experiments = polymer.get_number_of_experiments()
    logger.info('Loaded {} ({} Experiments)'.format(osp.join(path, sdf_file), nr_experiments))
    polymer.rephase_fids()
    return polymer

def create_plot_1_and_2(source_fid, source_mag_dec):
    logger.debug('creating plot 1 and 2')
    # create and plot figure 1
    p1 = figure(plot_width=800, plot_height=500,
                title='Free Induction Decay', webgl=True,
                lod_factor=1000,lod_interval=150,lod_threshold=1000)
    # convert data to handle in bokeh
    # ignore rephased_fid at this point
#     source_fid = ColumnDataSource(data=ColumnDataSource.from_df(rephased_fid))
    p1.line('index', 'im', source=source_fid, color='blue')
    p1.line('index', 'real', source=source_fid, color='green')
    p1.line('index', 'magnitude', source=source_fid, color='red')

    # create and plot figure 2
    p2 = figure(plot_width=400, plot_height=400,
                title='Magnetization Decay')
    p2.circle_cross('tau', 'phi_normalized', source=source_mag_dec, color="navy")
    p2.line('tau', 'fit_phi', source=source_mag_dec, color="teal")
    return p1, p2

def fit_mag_decay_all(polymer, par_df):
    logger.debug('fitting all mag decay')
    p3 = figure(plot_width=400, plot_height=400,
            title='normalized phi vs normalized tau', webgl = True,
                y_axis_type = 'log',
                x_axis_type = 'linear')

    nr_experiments = polymer.get_number_of_experiments()
    r1=np.zeros(nr_experiments)
    MANY_COLORS = 0
    p3_line_glyph=[]
    for i in range(1, nr_experiments):
        try:
            par=polymer.getparameter(i)
            fid=polymer.getfid(i)
            tau= polymer.get_tau_axis(i) # tau axis could be parameter created right after import
            try:
                startpoint=polymer.getparvalue(i,'fid_range')[0]
                endpoint=polymer.getparvalue(i,'fid_range')[1] 
            except:
                startpoint=int(0.05*polymer.getparvalue(i,'BS'))
                endpoint = int(0.1*polymer.getparvalue(i,'BS'))
                logger.warning('using preset range (from {startpoint} to {endpoint}) in experiment {i}'.format(**locals()))
                polymer.addparameter(i,'fid_range',(startpoint,endpoint))
            phi = get_mag_amplitude(fid, startpoint, endpoint,
                                    polymer.getparameter(i)['NBLK'], polymer.getparameter(i)['BS'])
            df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
            df['phi_normalized']=(df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] )
            polymer.addparameter(i,'df_magnetization',df)
            
            # initial values for the fit
            p0 = [1, 2 * polymer.getparvalue(i,'T1MX'), 0]
            df, popt = magnetization_fit(df,p0, fit_option=2)
            polymer.addparameter(i,'amp',popt[0])
            polymer.addparameter(i,'r1',popt[1])
            polymer.addparameter(i,'noise',popt[2])
            r1[i] = popt[1]
            tau = popt[1]*df.tau
            phi = popt[0]**-1*(df.phi_normalized - popt[2])
            p3_df=pd.DataFrame(data=np.c_[ tau, phi ], columns=['tau', 'phi'])
            source_p3=ColumnDataSource(data=ColumnDataSource.from_df(p3_df))
            p3_line_glyph.append(p3.line('tau', 'phi', source=source_p3))
            MANY_COLORS+=1
        except KeyError:
            logger.warning('no relaxation experiment found')
            polymer.addparameter(i,'amp',float('NaN'))
            polymer.addparameter(i,'r1',float('NaN'))
            polymer.addparameter(i,'noise',float('NaN'))
    COLORS=viridis(MANY_COLORS)
    for ic in range(MANY_COLORS):
        p3_line_glyph[ic].glyph.line_color=COLORS[ic]
    par_df['r1']=r1
    return p3

def modify_doc(doc):
    """ Contains the application, including all callbacks
        TODO: could the callbacks be outsourced?
    :param doc:
    :type doc:
    """
    logger.debug('modify_doc has been called')
    def load_polymer(filename='glyzerin_d3_300K.sdf', path=None):
        # load this datafile initially:
        polymer = load_data()
        start_ie = 1     # initially set ie = 1
        par_df, columns, discrete, continuous, time, quantileable = polymer.scan_parameters(20)
        # for the initial call get the dataframes without callback
        # they are being updated in following callbacks
        fid, df = get_data_frames(start_ie)
        source_fid = ColumnDataSource(data=ColumnDataSource.from_df(fid))
        source_mag_dec = ColumnDataSource(data=ColumnDataSource.from_df(df))
        return polymer, source_fid, source_mag_dec, par_df, columns, discrete, continuous, time, quantileable, start_ie

    def get_data_frames(ie,):
            """ Called one time initially, and then every time the experiment number is changed by the slider
            :param ie: experiment number
            :type ie: int
            :returns: dataframe from stella datafile and dataframe with tau and phi and fitted values
            :rtype: list of 2 pandas dataframes
            """
            logger.debug('get_dataframe with ie={}'.format(ie))
            fid = polymer.getfid(ie) #read FID or series of FIDs for selected experiment
            try:
                tau = polymer.get_tau_axis(ie) #numpy array containing the taus for experiment ie
                try:
                    startpoint=fid_slider.range[0] #lower integration bound
                    endpoint = fid_slider.range[1] #upper integration bound
                except NameError:
                    # fid_slider not initialized for first plot. Use default values:
                    startpoint=int(0.05*polymer.getparvalue(ie,'BS'))
                    endpoint = int(0.1*polymer.getparvalue(ie,'BS'))
                    logger.debug('fid_slider not initialized for first plot. Use default values {} and {}.'.format(startpoint, endpoint))
                    
                polymer.addparameter(ie,'fid_range',(startpoint,endpoint)) #add integration range to parameters to make it accesible
                phi = get_mag_amplitude(fid, startpoint, endpoint,
                                        polymer.getparvalue(ie,'NBLK'),
                                        polymer.getparvalue(ie,'BS')) # list containing averaged fid amplitudes (which is proportional to a magnetization phi)
                df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi']) # DataFrames are nice
                df['phi_normalized'] = (df['phi'] - df['phi'].iloc[0] ) / (df['phi'].iloc[-1] - df['phi'].iloc[1] ) #Normalize magnetization,
                #Note: in the normalized magnetization the magnetization build-up curves and magnetization decay curves look alike
                #Note: this makes it easier for fitting as everything looks like 1 * exp(-R/time) in first order
                polymer.addparameter(ie,'df_magnetization',df) # make the magnetization dataframe accesible as parameter
                fit_option = 2 #mono exponential, 3 parameter fit
                p0=[1.0, polymer.getparvalue(ie,'T1MX')**-1*2, 0] #choose startparameters for fitting an exponential decay
                df, popt = magnetization_fit(df, p0, fit_option) # use leastsq to find optimal parameters
                polymer.addparameter(ie,'popt(mono_exp)',popt) # add fitting parameters for later access
                logger.info('fitfunction(t) = {} * exp(- {} * t) + {}'.format(*popt)) # print the fitting parameters to console (for convenience)
            except KeyError:
                logger.warning('no relaxation experiment found')
                tau=np.zeros(1)
                phi=np.zeros(1)
                df = pd.DataFrame(data=np.c_[tau, phi], columns=['tau', 'phi'])
                df['phi_normalized'] = np.zeros(1)
                df['fit_phi'] = np.zeros(1)
            return fid, df

    def calculate_mag_dec(attr, old, new, start_ie=None):
        ''' Is being call from the callback for the experiment chooser
            loads selected experiment visualize in plot p1 and p2 
            gets experiment number from the slider
            writes source_fid.data from the fid from the polymer object
            writes source_mag_dec.data from the dataframe
            '''
        ie = experiment_slider.value   #get expermient number from the slider
        logger.debug('calculate mag_dec for ie={}'.format(ie))
        fid, df = get_data_frames(ie)
        source_fid.data=ColumnDataSource.from_df(fid) #convert fid to bokeh format
        source_mag_dec.data = ColumnDataSource.from_df(df)

    def plot_par():
        ''' Creates plot for the parameters 
            Called with every update from the callback'''
        logger.debug('creating plot for the parameters')

        # read data due to selection of select_x/y
        xs = par_df[select_xaxis.value ].values
        ys = par_df[select_yaxis.value].values
        # read titles due to name of select_x/y
        x_title = select_xaxis.value.title()
        y_title = select_yaxis.value.title()

        # remark: many attributes in a bokeh plot cannot be modified after initialization
        #         for example p4.x_axis_type='datetime' does not work. keywords are a
        #         workaround to pass all optional arguments initially
        # set optional keyword arguments, kw, for figure()
        kw = dict() #initialize
        if select_xaxis.value in discrete:
            kw['x_range'] = sorted(set(xs))
        if select_yaxis.value in discrete:
            kw['y_range'] = sorted(set(ys))
        if select_yaxis.value in time:
            kw['y_axis_type'] = 'datetime'
        if select_xaxis.value in time:
            kw['x_axis_type'] = 'datetime'
        kw['title']="%s vs %s" % (x_title, y_title)
        # create figure using optional keywords kw
        p4 = figure(plot_height=300, plot_width=600, tools='pan,box_zoom,reset',
                    **kw)
        # set axis label
        p4.xaxis.axis_label = x_title
        p4.yaxis.axis_label = y_title

        # strings at x axis ticks need a lot of space. solution: rotate label orientation
        if select_xaxis.value in discrete:
            p4.xaxis.major_label_orientation = pd.np.pi / 4 # rotates labels...

        # standard size of symbols
        sz = 9
        # custom size of symbols due to select_size
        if select_size.value != 'None':
            groups = pd.qcut(pd.to_numeric(par_df[select_size.value].values), len(SIZES))
            sz = [SIZES[xx] for xx in groups.codes]

        # standard color
        c = "#31AADE"        
        # custom color due to select_color
        if select_color.value != 'None':
            groups = pd.qcut(pd.to_numeric(par_df[select_color.value]).values, len(COLORS))
            c = [COLORS[xx] for xx in groups.codes]

        # create the plot using circles
        p4.circle(x=xs, y=ys, color=c, size=sz, line_color="white", alpha=0.6, hover_color='white', hover_alpha=0.5)
        return p4 #return the plot
    
    def callback_update_plot_1(attr, old, new):
        ''' Callback for update of figure 1 in parameters tab '''
        tabs.tabs[1].child.children[1] = plot_par()
        print(tabs.tabs[1].child.children[1])
        logger.debug('Parameter plot updated')
#        p4 = plot_par()

    def callback_update_p3():
        logger.debug('update plot 3')
        p3 = fit_mag_decay_all(polymer,par_df)
        return p3

    def callback_update_experiment(attr, old, new):
        """ Callback for the experiment chooser
        """
        ie = experiment_slider.value
        logger.debug('Callback experiment update, ie={}'.format(ie))
        fid_slider.end = polymer.getparvalue(ie,'BS')
        try:
            fid_slider.range=polymer.getparvalue(ie,'fid_range')
        except:
            startpoint = int(0.05 * polymer.getparvalue(ie,'BS'))
            endpoint = int(0.1 * polymer.getparvalue(ie,'BS'))
            fid_slider.range=(startpoint,endpoint)
        calculate_mag_dec(attr,old,new)
        
    def callback_load_more_data(attr,old,new):
        ''' callback for loading of data '''
        # first clear the list, or it keeps on growing with the same filenames:
        sdf_list = []
        path=pathbox.value.strip()
        file=filebox.value.strip()
        logger.info('callback for loading data. pathname: {}, filename: {}'.format(path, file))
        if file=="*.sdf":
            allsdf=filter(lambda x: x.endswith('.sdf'),os.listdir(path))
            for f in allsdf:
                sdf_list.append(sdf.StelarDataFile(f,path))
        else:
            sdf_list.append(sdf.StelarDataFile(file,path))
        
        filenames=[x.file() for x in sdf_list]
        filenames_df=pd.DataFrame(data=filenames,columns=['file'])
        table_source.data=ColumnDataSource.from_df(filenames_df)
        
    def callback_select_filename(attr,old,new):
        ''' callback for selecting one datafile '''
        logger.info('callback for selecting new filename: {}, {}, {}'.format(attr, old, new))

    def callback_export_data(attr,old,new):
        logger.debug('callback_export_data has been called ')
        logger.error('Not implemented!')
        pass
    
    def callback_write_table_to_file(attr,old,new): ##FIXME
        logger.debug('callback_write_table_to_file has been called ')
        logger.error('Not implemented!')
        pass
#        path=export_text.value.strip()
#        exportdata=export_source.data
#        CustomJS(args=dict(source=export_source),
#                 code=open(join(dirname(__file__), "export_csv.js")).read())

    def callback_update_parameters():
        ''' callback for button
            function to call when button is clicked
            for updates parameters of polymer, since they can change during evaluation '''
        logger.debug('callback for button (update parameter).')
        par_df, columns, discrete, continuous, time, quantileable = polymer.scan_parameters()
        select_xaxis.options=columns
        select_yaxis.options=columns
        select_size.options=['None']+quantileable
        select_color.options=['None']+quantileable

    logger.info('Starting the script')
    ### This is the start of the script ###
    ### The callbacks are above ###

    #load data:
    # TODO: how to handle multiple datafiles?
    # New Tab for each datafile?
    # dropdown selection to choose datafile
    # complete new start of process? (probably not prefered)
    polymer, source_fid, source_mag_dec, par_df, columns, discrete, continuous, time, quantileable, start_ie = load_polymer()
    
    # initialy creates the plots p1 and p2
    p1, p2 = create_plot_1_and_2(source_fid, source_mag_dec)
    
    ### initiates widgets, which will call the callback on change ###
    # initiate slider to choose experiment
    experiment_slider = Slider(start=1, end=polymer.get_number_of_experiments(), value=1, step=1,callback_policy='mouseup', width=800) #select experiment by value
    # initiate slider for the range in which fid shall be calculated
    # select the intervall from which magneitization is calculated from fid
    fid_slider = RangeSlider(start=1,end=polymer.getparvalue(start_ie,'BS'),
                             range=polymer.getparvalue(start_ie,'fid_range'),
                             step=1,callback_policy='mouseup', width=400)

    # fit magnetization decay for all experiments
    p3 = fit_mag_decay_all(polymer, par_df)
    # refit mag dec with updated ranges after button push
    button_refit = Button(label='Update',button_type="success")
    button_refit.on_click(callback_update_p3)

    # initialize empty source for experiment slider
    source = ColumnDataSource(data=dict(value=[]))
    # 'data' is the attribute. it's a field in source, which is a ColumnDataSource
    # initiate callback_update_experiment which is the callback
    source.on_change('data',callback_update_experiment) #source for experiment_slider
    experiment_slider.callback = CustomJS(args=dict(source=source),code="""
        source.data = { value: [cb_obj.value] }
    """)#unfortunately this customjs is needed to throttle the callback in current version of bokeh

    # initialize empty source for fid slider, same as above
    source2 = ColumnDataSource(data=dict(range=[], ie=[]))
    source2.on_change('data',calculate_mag_dec)
    fid_slider.callback=CustomJS(args=dict(source=source2),code="""
        source.data = { range: cb_obj.range }
    """)#unfortunately this customjs is needed to throttle the callback in current version of bokeh

    # same for the update button
    button_scan = Button(label='Scan Parameters',button_type="success")
    button_scan.on_click(callback_update_parameters)
    
    # here comes for callbacks for x, y, size, color
    select_xaxis = Select(title='X-Axis', value='ZONE', options=columns)
    select_xaxis.on_change('value', callback_update_plot_1)

    select_yaxis = Select(title='Y-Axis', value='TIME', options=columns)
    select_yaxis.on_change('value', callback_update_plot_1)

    select_size = Select(title='Size', value='None', options=['None'] + quantileable)
    select_size.on_change('value', callback_update_plot_1)

    select_color = Select(title='Color', value='None', options=['None'] + quantileable)
    select_color.on_change('value', callback_update_plot_1)

    controls_p4 = widgetbox([button_scan, select_xaxis,select_yaxis,select_color,select_size], width=150)
    #p4 = plot_par()
    layout_p4 = row(controls_p4,plot_par())
    logger.debug('layout for parameter plot created')

    ####
    #### TODO: write file input
    #### TODO: select files to import
    #### TODO: discard imported files
    ####

    # load more data:
    table_source=ColumnDataSource(data=dict())
    sdf_list=[polymer]
    filenames=[x.file() for x in sdf_list]
    files_df=pd.DataFrame(data=filenames,columns=['file'])
    table_source.data=ColumnDataSource.from_df(files_df)
    t_columns = [
        TableColumn(field='file', title='Path / Filename'),
        #TableColumn(field='file', title='Filename'),
        ]
    table=DataTable(source=table_source,columns=t_columns)
    pathbox=TextInput(title="Path",value=osp.join(osp.curdir, '/data'))
    filebox=TextInput(title="Filename",value="*.sdf")
    pathbox.on_change('value',callback_load_more_data)
    filebox.on_change('value',callback_load_more_data)
    layout_input=column(pathbox,filebox,table)
    # select a filename by clicking on it:
    table_source.on_change('selected', callback_select_filename)
    # Data Out: export data from figures
    #         & export parameters

    export_source=ColumnDataSource(data=dict())
    export_columns=[]
    output_table=DataTable(source=export_source,columns=export_columns)
    export_slider = Slider(start=1, end=4, value=3, step=1,callback_policy='mouseup', width=200) #do we need mouseup on this?
    export_slider.on_change('value',callback_export_data)
    export_text = TextInput(title="Path",value=os.path.curdir)
    export_button = Button(label='Export to csv',button_type="success") # FIXME Callback  doesn't work yet
    export_button.on_click(callback_write_table_to_file)
 
    layout_output=row(column(export_slider,export_text,export_button),output_table)
    logger.debug('after layout_output')
    

    # set the layout of the tabs
    layout_p1 = column(experiment_slider, p1,
                       row(
                           column(fid_slider,p2),
                           column(button_refit, p3)
                           ),
                       )
    tab_relaxation = Panel(child = layout_p1, title = 'Relaxation')
    tab_parameters = Panel(child = layout_p4, title = 'Parameters')
    tab_input = Panel(child = layout_input, title = 'Data In')
    tab_output = Panel(child = layout_output, title = 'Data Out')

    # initialize tabs object with 3 tabs
    tabs = Tabs(tabs = [tab_relaxation, tab_parameters,
                        tab_input, tab_output])
    logger.debug('tabs')
    doc.add_root(tabs)
    doc.add_root(source) # i need to add source to detect changes
    doc.add_root(source2)
    logger.debug('tab tab')

def main():
    bokeh_app = Application(FunctionHandler(modify_doc))
    io_loop = IOLoop.current()
    logger.info('Opening Bokeh application on http://localhost:5006/')
    server = Server({'/': bokeh_app}, io_loop=io_loop)
    server.start()
    io_loop.add_callback(server.show, "/")
    io_loop.start()

if __name__ == '__main__':
    main()

