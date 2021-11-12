

import pandas as pd
import matplotlib.pyplot as plt
import sys
import getopt
import math
import numpy as np
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import os



usage = """

*****************************************************************************
This script processes data obtained from pose-estimation of mouse behaviour
in an Open Maze test, using DeepLabCut.

Use the following command to run the analysis:
python processDLC.py <option-argument> <DLC-processed .h5 datafile>

options:
opt:			argument:	Default		description:
  -r (--resolution)	  string  	 1280x960 	  Resolution of video from which data was obtained
  -f (--framerate)	  integer  	 25 	  	  Framerate of video from which data was obtained
  -s (--scale)  	  float  	 6 	  	  Scale of the video.
  -t (--topleft)	  string  	 312x156 	  Coordinate of the top-left most corner of setup
  -b (--botright)	  string  	 978x839 	  Coordinate of the bottom-right most corner of setup
  -o (--outfile)	  string  	 ./DLCprocess.tsv Path to output filename
  -a (--append)           NA  	         False	  	  Whether to append DLC data
  -h (--help)		  NA  		 False	 	  Displays usage of the script
*****************************************************************************

"""


########################################
############ Class and Def #############
########################################

### function to parse command-line arguments
def argparse():
	# attempt to read command-line option-argument tuples and mandatory argument.
	try:
	    
		options, inputfile = getopt.getopt(sys.argv[1:], "hr:f:s:t:b:o:a",["help","framerate=", "scale=", "resolution=", "botright=", "topleft=", "outfile=", "append"])
		inputfile = str(inputfile[0]).strip()
		ext = os.path.splitext(inputfile)[1]
		if ext != ".h5":
		    sys.exit("Wrong input format provided")
	except getopt.GetoptError as err:
		sys.exit(usage)		# exit script and print usage if arguments to options are not provided
	except IndexError:
		sys.exit(usage)		# exit script and print usage if command-line input is not provided
 
    # define global variables
	global framerate
	global res
	global scale
	global topleft
	global botright
	global outfile
	global append
	
	# assign variables
	framerate = 25
	res = [1280,960]
	scale = 6
	topleft = [312,156]
	botright = [978,839]
	outfile="./DLCprocess.tsv"
	append = False
	
	# parse command-line options into its appropriate variables/actions
	for opt,arg in options:
		if opt in ("-r","--resolution"):
			res = list(map(int,arg.split('x')))
		if opt in ("-f","--framerate"):
			framerate = int(arg)
		if opt in ("-s","--scale"):
			scale = int(arg)
		if opt in ("-t","--topleft"):
			topleft = list(map(int,arg.split('x')))
		if opt in ("-b","--botright"):
			botright = list(map(int,arg.split('x')))
		if opt in ("-o","--outfile"):
			outfile = arg
		if opt in ("-a","--append"):
			append = True
		if opt in ("-h","--help"):
			sys.exit(usage)
	
	return inputfile

def runDLCprocess(path_to_h5, 
    frame_rate = 25, 
    image_scale = 6, 
    image_resolution = [1280,960], 
    topleft_pixel = [312,156], 
    botright_pixel = [978,839], 
    output_file = "./DLCprocess.tsv", 
    append_tsv = False):
    """
    Run DLCprocess workflow

    Parameters
    ----------
    path_to_h5 : str, required
        Path to input h5 file.
    frame_rate : integer, optional
        Frame rate of video from which data was obtained. The default is 25.
    image_scale : integer, optional
        Scale factor of setup relative to movie. The default is 6.
    image_resolution : 2-element integer array, optional
        Resolution of video. The default is [1280,960].
    topleft_pixel : 2-element integer array, optional
        x,y pixel position of the top-left position of the maze. The default is [312,156].
    botright_pixel : 2-element integer array, optional
        x,y pixel position of the bottom-right position of the maze. The default is [978,839].
    output_file : str, optional
        Path to output tsv file. By default, data will be exported to file named 'DLCprocess.tsv' in current working directory.
    append_tsv : bool, optional
        Whether or not to append information to output tsv. The default is false.

    Returns
    -------
    4 Matplotlib plots pertaining to locomotion, velocity, angular velocity and time spent in quadrants.
    A tsv file containing data for each mice. 

    """

    # define global variables
    global framerate
    global res
    global scale
    global topleft
    global botright
    global outfile
    global append

    # add info
    framerate = frame_rate
    res = image_resolution
    scale = image_scale
    topleft = topleft_pixel
    botright = botright_pixel
    outfile = output_file
    append = append_tsv

    runworkflow(path_to_h5)

    

    

### function to process input h5 file and extract info
def processh5(inputfile):
    # read data
    dat = pd.read_hdf(inputfile)
    
    # define global variables
    global parts
    global exp
    global mice
    global micepos
    
    # get experiment, mice names and body parts
    cols = dat.columns.to_frame()
    
    ## check column names
    desired_columns = ['scorer', 'individuals', 'bodyparts', 'coords']
    if desired_columns != list(cols.columns):
        sys.exit("Data do not contain correct headers")
    
    ## extract column info
    exp = list(set(cols['scorer']))[0]
    mice = list(set(cols['individuals']))
    mice.sort()
    parts = list(set(cols['bodyparts']))
    
    # exit if required parts are not labelled
    if not all(item in parts for item in ['middle_head','tail_base']):
        sys.exit("Data is missing coordinates for key bodyparts")
    
    # get mice positions in setup
    firstpos = np.array(dat.loc[0:0, (exp, mice, "nose", ["x","y"])])[0]
    firstpos = np.reshape(firstpos, (4,2))
    firstposdf = pd.DataFrame(firstpos, columns=['x','y'])
    firstposdf['mice'] = mice
    firstposdf['xdiff'] = (firstposdf['x'] - firstposdf['x'].mean()) > 0
    firstposdf['ydiff'] = (firstposdf['y'] - firstposdf['y'].mean()) > 0
    firstposdf['sum'] = firstposdf['xdiff'] | firstposdf['ydiff']
    
    for index, row in firstposdf.iterrows():
        if row['sum'] == False:
            firstposdf.loc[index, 'sum'] = "Top-left"
        elif row['xdiff'] == True & row['ydiff'] == True:
            firstposdf.loc[index, 'sum'] = "Bottom-right"
        elif row['xdiff'] == False and row['ydiff'] == True:
            firstposdf.loc[index, 'sum'] = "Bottom-left"
        elif row['xdiff'] == True and row['ydiff'] == False:
            firstposdf.loc[index, 'sum'] = "Top-right"
    micepos = list(firstposdf['sum'])
   
    # return datafile
    return dat

### function to test locomotion
def testloco(inputdat):
    
    # Get vector between the coord of middle_head between frames
    xdiff = np.diff(np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'middle_head', 'x')]), axis=0)
    ydiff = np.diff(np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'middle_head', 'y')]), axis=0)
    
    # Calculate displacement of vector and create output
    out = np.sqrt((xdiff **2) + (ydiff**2)) * 0.026458 * scale * 0.01
    out = np.r_[np.zeros((1,4)), out]
    out = pd.DataFrame(out, columns = micepos)

    return(out)
    
### function to test velocity
def testvelo(dat):
    vel = pd.DataFrame(np.array(dat)*25, columns = micepos)
    vellroll = vel.rolling(25).mean().shift(-3)
    vellroll.fillna(0)

    return(vellroll)

# function to test rotation
def testrot(inputdat):
    
    # extract vector of middle-head to tail 
    ax = np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'tail_base', 'x')])
    ay = np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'tail_base', 'y')])
    bx = list(np.array(inputdat.loc[0:0, (exp, mice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:0, (exp, mice, 'tail_base', 'x')])) + list(np.array(inputdat.loc[0:(len(inputdat)-2), (exp, mice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:(len(inputdat)-2), (exp, mice, 'tail_base', 'x')]))
    by = list(np.array(inputdat.loc[0:0, (exp, mice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:0, (exp, mice, 'tail_base', 'y')])) + list(np.array(inputdat.loc[0:(len(inputdat)-2), (exp, mice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:(len(inputdat)-2), (exp, mice, 'tail_base', 'y')]))

    # count angle change
    outraw = np.arctan2(np.array(ay), np.array(ax)) - np.arctan2(np.array(by), np.array(bx))
    
    # normalize values so that 2pi is represented as 1 rotation
    ## Create vectors with 2pi
    twopie = np.full((len(inputdat), len(mice)), (math.pi*2))
    
    ## normalize
    outraw = np.where(outraw > math.pi, outraw - (2*math.pi), outraw)
    outraw = np.where(outraw < -math.pi, outraw + (2*math.pi), outraw)
    outcumsum = np.cumsum(outraw, axis=0)
    outcumsum = np.divide(outcumsum, twopie)
    
    # Count number of rotations
    ## make pd dataframe for cumsum data and new outout dataframe
    out = pd.DataFrame(outcumsum, columns = micepos)
    norm = pd.DataFrame(np.zeros((len(inputdat), len(mice))), columns = micepos)
    for thismice in micepos:
        thresh = 0
        for index, row in out.iterrows():
            if row[thismice] > (thresh+1):
                thresh+=1
            elif row[thismice] < (thresh-1):
                thresh-=1
            norm.loc[index, thismice] = row[thismice] - thresh
            
    # get change in angle between frames and select timepoint when change is above 0.9
    outdiff = norm.diff()
    outdiff['time'] = outdiff.index/framerate
    outdiffroll = pd.melt(outdiff, value_vars=micepos, id_vars='time')
    outdiffsigroll = outdiffroll[(outdiffroll['value'] > 0.9) | (outdiffroll['value'] < -0.9)]
    
    # Calculate angular velocity
    outrawdf = pd.DataFrame((outraw*25)/(math.pi), columns = micepos)
    angveldf = outrawdf.rolling(25).mean().shift(-3)
    angveldf.fillna(0)

    return(outdiffsigroll, angveldf)

    

def testquad(inputdat):
    
    # get coordiantes of boundaries
    xbounds = list(map(lambda x: (topleft[0] + ((botright[0]-topleft[0]-25)/8)*(x+1)) if x < 4 else (topleft[0] + 25 + ((botright[0]-topleft[0]-25)/8)*(x+1)), range(8)))
    ybounds = list(map(lambda x: (topleft[1] + ((botright[1]-topleft[1]-25)/8)*(x+1)) if x < 4 else (topleft[1] + 25 + ((botright[1]-topleft[1]-25)/8)*(x+1)), range(8)))  
  
    
    # x1 = topleft[0] + ((botright[0]-topleft[0]-25)/8)*1
    # x1 = topleft[0] + (botright[0]-topleft[0]-25)/
    # xmid = topleft[0] + (botright[0]-topleft[0])/2
    # x2 = botright[0] - (botright[0]-topleft[0]-25)/4
    # 
    # y1 = topleft[1] + (botright[1]-topleft[1]-25)/4
    # ymid = topleft[1] + (botright[1]-topleft[1])/2
    # y2 = botright[1] - (botright[1]-topleft[1]-25)/4
    
    
    # prepare pos of basetail
    basexcoord = np.array(inputdat.loc[0:len(inputdat), (exp, mice, "tail_base", ["x"])])
    baseycoord = np.array(inputdat.loc[0:len(inputdat), (exp, mice, "tail_base", ["y"])])
    
    # test x axis
    xres = np.zeros((len(inputdat), len(mice)))
    for bound in xbounds:
        xres = xres + (basexcoord > bound) 
        
    # test y axis
    yres = np.zeros((len(inputdat), len(mice)))
    for bound in ybounds:
        yres = yres + (baseycoord > bound)    
        
    xres = xres + 1
    yres = yres * 8
    totalres = xres + yres
    # totalres = np.char.mod('%d', totalres)
    
    totalrespos = np.full((len(inputdat), len(mice)), 'None')
    corners = [1,5,33,37,4,8,36,40,25,29,57,61,28,29,60,64]
    middle = [10,11,18,19,14,15,22,23,42,43,46,47,50,51,54,55]
    totalrespos = np.where(np.isin(totalres, corners), "Corner", totalrespos)
    totalrespos = np.where(np.isin(totalres, middle), "Middle", totalrespos)
    
    
    # totalrespos = np.where(np.isin(totalres, [25,29,57,61]), "Bottom-left", totalrespos)
    # totalrespos = np.where(np.isin(totalres, [28,29,60,64]), "Bottom-right", totalrespos)
    
    
    # tally data
    totalresposdf = pd.DataFrame(totalrespos, columns = micepos)
    totalresposdf['id'] = totalresposdf.index
    totalresposdfmelt = pd.melt(totalresposdf, value_vars=micepos, id_vars='id')
    totalresposdfmelt = totalresposdfmelt[totalresposdfmelt["value"] != "None"]
    totalresposdfmelt['size'] = totalresposdfmelt.groupby(['variable','value']).transform(np.size)
    
    # remove duplicates and prepare output dataframe
    totalresposdfmeltunique = totalresposdfmelt.loc[0:totalresposdfmelt.index[-1], ('variable', 'value','size')].drop_duplicates()
    totalresposdfmeltunique = totalresposdfmeltunique.sort_values(by = ['variable','value'])
    
    totalresposdfmeltunique['cumsum'] = totalresposdfmeltunique.groupby(['variable']).cumsum()  # this convert counts to cumulative sum
    
    out = totalresposdfmeltunique.pivot_table(index='variable',columns='value',values='cumsum', fill_value = 0)

    newdf = pd.DataFrame(np.array(out)/len(inputdat), columns = list(out.columns))
    newdf['mice'] = out.index
    return(newdf)

def testrear(inputdat):
    
     # extract vector of middle-head to tail 
    ax = np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'tail_base', 'x')])
    ay = np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'tail_base', 'y')])
    
    axnose = np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'nose', 'x')]) - np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'tail_base', 'x')])
    aynose = np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'nose', 'y')]) - np.array(inputdat.loc[0:len(inputdat), (exp, mice, 'tail_base', 'y')])
    
    # Calculate displacement of vector and create output
    out = np.sqrt((ax **2) + (ay**2)) * 0.026458 * scale * 0.01
    out = np.r_[np.zeros((1,4)), out]
    out = pd.DataFrame(out, columns = micepos)
    out = out.rolling(25).mean().shift(-3)
    out.fillna(0)
    
    out2 = np.sqrt((axnose **2) + (aynose **2)) * 0.026458 * scale * 0.01
    out2 = np.r_[np.zeros((1,4)), out2]
    out2 = pd.DataFrame(out2, columns = micepos)
    out2 = out2.rolling(25).mean().shift(-3)
    out2.fillna(0)
    
    out2 = out2-out

    return(out,out2)

def showplot(loco, vellroll, angvel, quad, rear,rear2):
    newmicepos = ['Top-left', 'Top-right', 'Bottom-left', 'Bottom-right']

    # prepare and plot displacement data
    lococumsum = loco.cumsum()
    lococumsum['time'] = loco.index/framerate
    lococumsum = pd.melt(lococumsum, value_vars=newmicepos, id_vars='time')
    lococumsum.columns = ['Time [s]', 'Mice', 'Distance travelled [m]']
    
    
    #fig = plt.figure(num='Plot of mice locomotion')
    fig = plt.figure()
    sns.lineplot(x="Time [s]", y="Distance travelled [m]",
             hue="Mice",
             data=lococumsum, palette = "hls")

    
    # prepare and plot velocity data
    
    vellroll = vellroll.reindex(columns = newmicepos)
    
    vellroll['time'] = vellroll.index/framerate
    vellroll = pd.melt(vellroll, value_vars=newmicepos, id_vars='time')
    vellroll.columns = ['Time [s]', 'Mice', 'Velocty [m/s]']
    
    # fig = plt.figure(num='Velocity')
    # sns.lineplot(x="Time [s]", y="Velocty [cm/s]",
    #          hue="Mice",
    #          data=vellroll)
    g = sns.FacetGrid(vellroll, col="Mice", palette = "hls", col_wrap = 2, hue = "Mice")
    g.map(sns.lineplot, 'Time [s]', 'Velocty [m/s]', ci = [100000]*len(vellroll))

    
    # prepare and plot angular velocity data
    angvel['time'] = angvel.index/framerate
    angvel = pd.melt(angvel, value_vars=newmicepos, id_vars='time')
    angvel.columns = ['Time [s]', 'Mice', 'Angular velocty [rad/s]']
    # fig = plt.figure('')
    # sns.lineplot(x="Time [s]", y="Angular velocty [rad/s]",
    #          hue="Mice", 
    #          data=angvel)
    vel = sns.FacetGrid(angvel, col="Mice", palette = "hls", col_wrap = 2, hue = "Mice")
    vel.map(sns.lineplot, 'Time [s]', 'Angular velocty [rad/s]', ci = [100000]*len(vellroll))
    
    # prepare and plot stacked data
    pal = sns.color_palette("hls")
    #fig = plt.figure(num='Plot of time spent in quadrants')
    fig = plt.figure()
    bar1 =  sns.barplot(x="mice",  y="Middle", data=quad, color=pal[0])
    bar2 =  sns.barplot(x="mice",  y="Corner", data=quad, color=pal[1])
    # bar2 =  sns.barplot(x="mice",  y="Bottom-right", data=quad, color=pal[2])
    # bar1 =  sns.barplot(x="mice",  y="Bottom-left", data=quad, color=pal[3])
    # add legend
    first_bar = mpatches.Patch(color=pal[0], label='Time spent in middle')
    second_bar = mpatches.Patch(color=pal[1], label='Time spent in corners')
    # third_bar = mpatches.Patch(color=pal[2], label='Bottom-right quadrant')
    # fourth_bar = mpatches.Patch(color=pal[3], label='Bottom-left quadrant')
    plt.legend(handles=[first_bar, second_bar], prop={'size':5})
    bar1.set_ylabel("Proportion spent in quadrants")
    bar1.set_xlabel("Mice")
    
    # # prepare and plot rearing data
    # rear = rear.reindex(columns = newmicepos)

    # rear['time'] = rear.index/framerate
    # rear = pd.melt(rear, value_vars=newmicepos, id_vars='time')
    # rear.columns = ['Time [s]', 'Mice', 'Horizontal length']

    # rearplot = sns.FacetGrid(rear, col="Mice", palette = "hls", col_wrap = 2, hue = "Mice")
    # rearplot.map(sns.lineplot, 'Time [s]', 'Horizontal length', ci = [100000]*len(rear))
    
    # rear2 = rear2.reindex(columns = newmicepos)

    # rear2['time'] = rear2.index/framerate
    # rear2 = pd.melt(rear2, value_vars=newmicepos, id_vars='time')
    # rear2.columns = ['Time [s]', 'Mice', 'Horizontal length']

    # rearplot2 = sns.FacetGrid(rear2, col="Mice", palette = "hls", col_wrap = 2, hue = "Mice")
    # rearplot2.map(sns.lineplot, 'Time [s]', 'Horizontal length', ci = [100000]*len(rear2))
    
    
    plt.show()
    

def createtab(loco, vellroll, angrot, angvel, quad):
    # prepare output dataframe
    out = pd.DataFrame()
    out['Experiment'] = [exp]*4
    out['Mice position'] = micepos
    
    # displacement
    out['Displacement [m]'] = list(loco.sum())
    
    # speed
    out['Mean velocity [m/s]'] = list(vellroll.mean())
    out['Max velocity [m/s]'] = list(vellroll.max())
    
    # angular
    # angrotclock = angrot[angrot['value'] < 0]
    # print(angrotclock)
    nrot = list(angrot.groupby(['variable']).size())
    nrot = pd.DataFrame(nrot)
    nrot.index =  ['Bottom-left', 'Bottom-right', 'Top-left', 'Top-right']
    nrot = nrot.loc[['Bottom-right', 'Top-left', 'Top-right', 'Bottom-left']]
    out['Num of rotations'] = list(nrot[0])
    
    out['Mean angular velocity [m/s]'] = list(abs(angvel).mean())
    out['Max angular velocity [m/s]'] = list(abs(angvel).max())
    
    
    # quadrants (in progoress)
    #newquad = quad.set_index('mice').T.diff()
    
    # append data if requested
    if append:
        # check if file exists
        if os.path.isfile(outfile):
            dat = pd.read_csv(outfile, sep='\t')
            out = dat.append(out, ignore_index=True)
        else:
            print("File to append `%s` not found, creating a new file" % (outfile))
    
    
    # export data
    path = outfile.split("/")
    path = list(filter(None, path))
    path.remove(".")
    ## test if absolute path is given
    if os.path.exists("/"+ path[0]):
        out.to_csv(outfile,index=False, sep = "\t")
    else:
        newpath = ["."] + path[:-1]
        creatpath = "/".join(newpath)
        os.makedirs(creatpath, exist_ok = True)
        out.to_csv(outfile,index=False, sep = "\t")
        
def runworkflow(pathtoinput):
    # parse command-line arguments and process input file 
    datatable = processh5(pathtoinput)
    
    # run function to quantify measurements
    locodat = testloco(datatable)
    velodat = testvelo(locodat)
    nrotations,angveldat = testrot(datatable)
    quaddat = testquad(datatable)
    reardat,reardat2 = testrear(datatable)
    
    # output tables and plots
    createtab(locodat, velodat, nrotations, angveldat,quaddat)
    showplot(locodat, velodat, angveldat,quaddat,reardat,reardat2)

    
    

########################################
################ Main ##################
########################################

#if __name__ == "__main__":
    
    
    








