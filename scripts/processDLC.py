

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
    x1 = topleft[0] + (botright[0]-topleft[0]-25)/4
    xmid = topleft[0] + (botright[0]-topleft[0])/2
    x2 = botright[0] - (botright[0]-topleft[0]-25)/4
    y1 = topleft[1] + (botright[1]-topleft[1]-25)/4
    ymid = topleft[1] + (botright[1]-topleft[1])/2
    y2 = botright[1] - (botright[1]-topleft[1]-25)/4
    
    
    # prepare pos of basetail
    basexcoord = np.array(inputdat.loc[0:len(inputdat), (exp, mice, "tail_base", ["x"])])
    baseycoord = np.array(inputdat.loc[0:len(inputdat), (exp, mice, "tail_base", ["y"])])
    
    # test x axis
    xres = np.zeros((len(inputdat), len(mice)))
    for bound in [x1,xmid,x2]:
        xres = xres + (basexcoord > bound) 
        
    # test y axis
    yres = np.zeros((len(inputdat), len(mice)))
    for bound in [y1,ymid,y2]:
        yres = yres + (baseycoord > bound)    
        
    xres = xres + 1
    yres = yres * 4
    totalres = xres + yres
    # totalres = np.char.mod('%d', totalres)
    
    totalrespos = np.full((len(inputdat), len(mice)), 'None')
    totalrespos = np.where(np.isin(totalres, [1,3,9,11]), "Top-left", totalrespos)
    totalrespos = np.where(np.isin(totalres, [2,4,10,12]), "Top-right", totalrespos)
    totalrespos = np.where(np.isin(totalres, [5,7,13,15]), "Bottom-left", totalrespos)
    totalrespos = np.where(np.isin(totalres, [6,8,14,16]), "Bottom-right", totalrespos)
    
    # tally data
    totalresposdf = pd.DataFrame(totalrespos, columns = micepos)
    totalresposdf['id'] = totalresposdf.index
    totalresposdfmelt = pd.melt(totalresposdf, value_vars=micepos, id_vars='id')
    totalresposdfmelt['size'] = totalresposdfmelt.groupby(['variable','value']).transform(np.size)
    
    # remove duplicates and prepare output dataframe
    totalresposdfmeltunique = totalresposdfmelt.loc[0:len(totalresposdfmelt), ('variable', 'value','size')].drop_duplicates()
    totalresposdfmeltunique = totalresposdfmeltunique.sort_values(by = ['variable','value'])
    totalresposdfmeltunique['cumsum'] = totalresposdfmeltunique.groupby(['variable']).cumsum()  # this convert counts to cumulative sum
    out = totalresposdfmeltunique.pivot_table(index='variable',columns='value',values='cumsum', fill_value = 0)
    
    newdf = pd.DataFrame(np.array(out)/len(inputdat), columns = list(out.columns))
    newdf['mice'] = out.index
    return(newdf)

  
    

def showplot(loco, vellroll, angvel, quad):

    # prepare and plot displacement data
    lococumsum = loco.cumsum()
    lococumsum['time'] = loco.index/framerate
    lococumsum = pd.melt(lococumsum, value_vars=micepos, id_vars='time')
    lococumsum.columns = ['Time [s]', 'Mice', 'Distance travelled [m]']
    
    fig = plt.figure(num='Locomotion')
    sns.lineplot(x="Time [s]", y="Distance travelled [m]",
             hue="Mice",
             data=lococumsum, palette = "hls")

    
    # prepare and plot velocity data
    newmicepos = ['Top-left', 'Top-right', 'Bottom-left', 'Bottom-right']
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
    fig = plt.figure(num='Quadrants')
    bar4 =  sns.barplot(x="mice",  y="Top-right", data=quad, color=pal[0])
    bar3 =  sns.barplot(x="mice",  y="Top-left", data=quad, color=pal[1])
    bar2 =  sns.barplot(x="mice",  y="Bottom-right", data=quad, color=pal[2])
    bar1 =  sns.barplot(x="mice",  y="Bottom-left", data=quad, color=pal[3])
    # add legend
    first_bar = mpatches.Patch(color=pal[0], label='Top-right quadrant')
    second_bar = mpatches.Patch(color=pal[1], label='Top-left quadrant')
    third_bar = mpatches.Patch(color=pal[2], label='Bottom-right quadrant')
    fourth_bar = mpatches.Patch(color=pal[3], label='Bottom-left quadrant')
    plt.legend(handles=[first_bar, second_bar, third_bar, fourth_bar], prop={'size':5})
    bar1.set_ylabel("Proportion spent in quadrants")
    bar1.set_xlabel("Mice")
    
    
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
        

    
    

########################################
################ Main ##################
########################################

if __name__ == "__main__":
    
    # parse command-line arguments and process input file 
    pathtoinput = argparse()
    datatable = processh5(pathtoinput)
    
    # run function to quantify measurements
    locodat = testloco(datatable)
    velodat = testvelo(locodat)
    nrotations,angveldat = testrot(datatable)
    quaddat = testquad(datatable)
    
    # output tables and plots
    createtab(locodat, velodat, nrotations, angveldat,quaddat)
    showplot(locodat, velodat, angveldat,quaddat)
    








