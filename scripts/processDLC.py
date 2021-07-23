

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
This script.....

Use Syntax
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

def argparse():
	# attempt to read command-line option-argument tuples and mandatory argument.
	try:
	    
		options, inputfile = getopt.getopt(sys.argv[1:], "hr:f:s:t:b:o:a",["help","framerate=", "scale=", "resolution=", "botright=", "topleft=", "outfile=", "append"])
		inputfile = str(inputfile[0]).strip()
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
	
	# define variables
	framerate = 25
	res = [1280,960]
	scale = 6
	topleft = [312,156]
	botright = [978,839]
	outfile="./DLCprocess.tsv"
	append = False
	
	

	# # parse command-line options into its appropriate variables/actions
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
	# 
	# # Additional preparations and checks
	# dict_codonrank,aa = pickle.load(open('definitions','r'))		#pre-load dictionaries and aa list
	# dict_species = pickle.load(open('species','r'))
	# if not re.search(r'.txt',inputfile):
	# 	if not set(inputfile.upper()) <= set(aa):sys.exit('Error: Input sequence may contain illegal characters')			# check for non-amino acid characters
	# if len(weights) != len(list_organisms): sys.exit('Error: weights do not match number of organisms')			# check for correct number of weights
	# 
	
	return inputfile


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
    exp = list(set(cols['scorer']))[0]
    mice = list(set(cols['individuals']))
    mice.sort()
    parts = list(set(cols['bodyparts']))
    
    # get mice positions
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
   
    
    
    # firstposdf = firstposdf.sort_values(by=['x'], ignore_index=True)
    # firstposdf['xpos'] = firstposdf.index
    # firstposdf = firstposdf.sort_values(by=['y'], ignore_index=True)
    # firstposdf['ypos'] = firstposdf.index
    # firstposdf = firstposdf.sort_values(by=['xpos','ypos'], ignore_index=True)
    
    
    # # convert dat to longformat
    # newcolnames = np.core.defchararray.add(list(cols['individuals']), (["-"] * len(cols.index)))
    # newcolnames = np.core.defchararray.add(newcolnames, list(cols['bodyparts']))
    # newcolnames = np.core.defchararray.add(newcolnames,  (["-"] * len(cols.index)))
    # newcolnames = np.core.defchararray.add(newcolnames, list(cols['coords']))
    # 
    # dat.columns = newcolnames
    # dat["frame"] = dat.index
    # datlong = pd.melt(dat, value_vars=newcolnames, id_vars=['frame'])
    # datlong['individual'] = list(val[0] for val in datlong['variable'].str.split(pat = "-"))
    # datlong['bodyparts'] = list(val[1] for val in datlong['variable'].str.split(pat = "-"))
    # datlong['coords'] = list(val[2] for val in datlong['variable'].str.split(pat = "-"))
    # datlong['id'] = np.core.defchararray.add(list(datlong['bodyparts']), list(datlong['coords']))
    # datfinal = datlong.pivot(index='id', columns = 'individual', values = 'value')
    # print(datfinal)
    
    return dat


def testloco(inputdat):
    xdf = pd.DataFrame(columns=mice)
    ydf = pd.DataFrame(columns=mice)
    prevxdf = pd.DataFrame(columns=mice)
    prevydf = pd.DataFrame(columns=mice)
    for thismice in mice:
        xdf[thismice] = list(inputdat.loc[0:len(inputdat), (exp, thismice, 'middle_head', 'x')])
        ydf[thismice] = list(inputdat.loc[0:len(inputdat), (exp, thismice, 'middle_head', 'y')])

        prevxdf[thismice] = list(inputdat.loc[0:0, (exp, thismice, 'middle_head', 'x')]) + list(inputdat.loc[0:(len(inputdat)-2), (exp, thismice, 'middle_head', 'x')])
        prevydf[thismice] = list(inputdat.loc[0:0, (exp, thismice, 'middle_head', 'y')]) + list(inputdat.loc[0:(len(inputdat)-2), (exp, thismice, 'middle_head', 'y')])

    out = np.sqrt(((np.array(xdf) - np.array(prevxdf)) **2) + ((np.array(ydf) - np.array(prevydf))**2)) * 0.026458 * scale * 0.01
    out = pd.DataFrame(out, columns = micepos)



    
    # 
    # # create output array
    # temp_array = np.empty((len(inputdat.index),4))
    # out = pd.DataFrame(columns=mice)
    # 
    # #create dict with starting head coord
    # dict_currentcoord= {}
    # for thismice in mice:
    #     dict_currentcoord[thismice] = [float(inputdat.loc[0,(exp, thismice, 'middle_head', 'x')]), float(inputdat.loc[0,(exp, thismice, 'middle_head', 'y')])]
    # 
    # for index, row in inputdat.iterrows():
    #     for miceindex in range(len(mice)):
    #         thismice = mice[miceindex]
    #         #get x and y cooord
    #         x = float(inputdat.loc[index, (exp, thismice, 'middle_head', 'x')])
    #         y = float(inputdat.loc[index, (exp, thismice, 'middle_head', 'y')])
    # 
    #         # get distance travelled
    #         dist = math.sqrt((x-dict_currentcoord[thismice][0])**2 + (y-dict_currentcoord[thismice][1])**2) * 0.26458 * scale
    #         temp_array[index][miceindex] = dist
    # 
    #         # update dict
    #         dict_currentcoord[thismice] = [float(inputdat.loc[index,(exp, thismice, 'middle_head', 'x')]), float(inputdat.loc[index,(exp, thismice, 'middle_head', 'y')])]

    #
    # for key in dict_currentcoord:
    #     print(key, ' : ', dict_currentcoord[key])

    # df = pd.DataFrame(temp_array, columns = mice)
    #print(df.head())
    #print(xdf)
    #print(ydf)
    return(out)
    

def testvelo(dat):
    vel = pd.DataFrame(np.array(dat)*25, columns = micepos)
    vellroll = vel.rolling(25).mean().shift(-3)
    vellroll.fillna(0)
    return(vellroll)

def testrot(inputdat):
    ax = pd.DataFrame(columns=mice)
    ay = pd.DataFrame(columns=mice)
    bx = pd.DataFrame(columns=mice)
    by = pd.DataFrame(columns=mice)
    for thismice in mice:
        ax[thismice] = np.array(inputdat.loc[0:len(inputdat), (exp, thismice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:len(inputdat), (exp, thismice, 'tail_base', 'x')])
        ay[thismice] = np.array(inputdat.loc[0:len(inputdat), (exp, thismice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:len(inputdat), (exp, thismice, 'tail_base', 'y')])
        
        bx[thismice] = list(np.array(inputdat.loc[0:0, (exp, thismice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:0, (exp, thismice, 'tail_base', 'x')])) + list(np.array(inputdat.loc[0:(len(inputdat)-2), (exp, thismice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:(len(inputdat)-2), (exp, thismice, 'tail_base', 'x')]))
        by[thismice] = list(np.array(inputdat.loc[0:0, (exp, thismice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:0, (exp, thismice, 'tail_base', 'y')])) + list(np.array(inputdat.loc[0:(len(inputdat)-2), (exp, thismice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:(len(inputdat)-2), (exp, thismice, 'tail_base', 'y')]))

        #bx[thismice] = np.repeat(np.array(inputdat.loc[0:0, (exp, thismice, 'middle_head', 'x')]) - np.array(inputdat.loc[0:0, (exp, thismice, 'tail_base', 'x')]), len(inputdat))
        #by[thismice] = np.repeat(np.array(inputdat.loc[0:0, (exp, thismice, 'middle_head', 'y')]) - np.array(inputdat.loc[0:0, (exp, thismice, 'tail_base', 'y')]), len(inputdat))

    outraw = np.arctan2(np.array(ay), np.array(ax)) - np.arctan2(np.array(by), np.array(bx))
    outraw = np.where(outraw > math.pi, outraw - (2*math.pi), outraw)
    outraw = np.where(outraw < -math.pi, outraw + (2*math.pi), outraw)
    outcumsum = np.cumsum(outraw, axis=0)
    twopie = np.full((len(inputdat), len(mice)), (math.pi*2))
    onemod = np.full((len(inputdat), len(mice)), 1)

    outcumsum = np.divide(outcumsum, twopie)
    # outcumsum = np.where(outcumsum > 0, outcumsum - np.floor(outcumsum), outcumsum)
    # outcumsum = np.where(outcumsum < 0, outcumsum - np.ceil(outcumsum), outcumsum)

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
    #print(norm.head())

    
    outdiff = norm.diff()
    outdiff['time'] = outdiff.index/framerate
    outdiffroll = pd.melt(outdiff, value_vars=micepos, id_vars='time')
    outdiffsigroll = outdiffroll[(outdiffroll['value'] > 0.9) | (outdiffroll['value'] < -0.9)]
    
    
    #peak_indices = find_peaks((np.array(out['individual4'])),height = 0.95)
    # idx = (np.argwhere(np.diff(np.sign(np.array([0]*len(inputdat)) - np.array(out['individual3'])))).flatten())/25
    # print(idx )
    outrawdf = pd.DataFrame((outraw*25)/(math.pi), columns = micepos)
    angveldf = outrawdf.rolling(25).mean().shift(-3)
    angveldf.fillna(0)
    # outroll = norm
    #out = out.cumsum()
    #outroll = out.rolling(25).mean().shift(-3)
    return(outdiffsigroll, angveldf)
    # prepare and plot displacement data
    # angveldf['time'] = angveldf.index/framerate
    # angveldf = pd.melt(angveldf, value_vars=mice, id_vars='time')
    # angveldf.columns = ['Time [s]', 'Mice', 'Angular velocity [rad/s]']
    # 
    # 
    # g = sns.FacetGrid(angveldf, col="Mice", palette = "Set1", col_wrap=2)
    # g.map(sns.lineplot, 'Time [s]', 'Angular velocity [rad/s]')
    # # sns.lineplot(x="Time [s]", y="Angle relative to time 0",
    # #          hue="Mice", 
    # #          data=outroll)
    # plt.show()
    

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
    
    totalresposdfmeltunique = totalresposdfmelt.loc[0:len(totalresposdfmelt), ('variable', 'value','size')].drop_duplicates()
    totalresposdfmeltunique = totalresposdfmeltunique.sort_values(by = ['variable','value'])
    totalresposdfmeltunique['cumsum'] = totalresposdfmeltunique.groupby(['variable']).cumsum()
    a = totalresposdfmeltunique.pivot_table(index='variable',columns='value',values='cumsum', fill_value = 0)
    
    newdf = pd.DataFrame(np.array(a)/len(inputdat), columns = list(a.columns))
    newdf['mice'] = a.index
    return(newdf)

    


    # print(totalresposdf.groupby(['variable','value']).transform(np.size))

    
    
    

def showplot(loco, vellroll, angvel, quad):
    # fig=plt.figure(figsize = (10,30))
    # gs=GridSpec(3,2)
    # ax1=fig.add_subplot(gs[0,:])
    # ax2=fig.add_subplot(gs[1,:])
    # ax3=fig.add_subplot(gs[2,:])
    #fig, axs = plt.subplots(nrows=2, figsize=(10,20))
    
    # prepare and plot displacement data
    lococumsum = loco.cumsum()
    lococumsum['time'] = loco.index/framerate
    lococumsum = pd.melt(lococumsum, value_vars=micepos, id_vars='time')
    lococumsum.columns = ['Time [s]', 'Mice', 'Distance travelled [m]']
    
    # 
    fig = plt.figure(num='Locomotion')
    sns.lineplot(x="Time [s]", y="Distance travelled [m]",
             hue="Mice",
             data=lococumsum, palette = "hls")
    # ax = datcumsum.plot()
    # ax.set_ylabel('Distance travelled [cm]')
    # ax.set_xlabel('Time [s]')
    
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
    
    
    # quadrants
    #newquad = quad.set_index('mice').T.diff()
    
    
    # export
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
    pathtoinput = argparse()
    datatable = processh5(pathtoinput)
    
    # run function to calculate locomotion
    locodat = testloco(datatable)
    velodat = testvelo(locodat)
    nrotations,angveldat = testrot(datatable)
    quaddat = testquad(datatable)
    createtab(locodat, velodat, nrotations, angveldat,quaddat)
    showplot(locodat, velodat, angveldat,quaddat)
    








