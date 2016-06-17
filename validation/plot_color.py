import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def color_outlines(bp, column_name, color):
    """Color outlines (of the box, median lines, whiskers and caps) with color
    
       @param bp             handle of boxplot
       @param column_name    column name of dataframe that contained boxplot data
       @param color          color to color outlines with
    """

    parts = ['boxes', 'whiskers', 'medians', 'caps', 'fliers']
    for part in parts:
        plt.setp(bp[column_name][part], color=color)
        if part=='fliers':
            for flier in bp[column_name]['fliers']:
                flier.set(marker='x', color=color, alpha=0.5)


def fill_boxes_with_color(numBoxes, bp, ax, column_name, color):
    """Fill boxplot boxes with a light shaded (alpha=0.5) color
    
       @param numBoxes       total number of boxes in boxplot (with handle bp on axes ax)
       @param bp             handle of boxplot
       @param ax             axes handle 
       @param column_name    column name of dataframe that contained boxplot data
       @param color          color to fill in boxes with
    """
    
    # color outlines first
    color_outlines(bp, column_name, color)


    from matplotlib.patches import Polygon
    medians = list(range(numBoxes))
    for i in range(numBoxes):
    
        # box in question
        box = bp[column_name]['boxes'][i]
        
        # get x,y position dat of box
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
     
        # use Polygon to fill in with 'color' with shading alpha=0.5
        boxPolygon = Polygon(boxCoords, facecolor=color, alpha=0.5)
        ax.add_patch(boxPolygon)
    
        # Now draw the median lines back over what we just filled in
        med = bp[column_name]['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            plt.plot(medianX, medianY, color)
            medians[i] = medianY[0]




def color_boxplot(ax, df1, df2, color_columns, redshift_columns, zmin, zmax, dbin, \
                  whis=1.5, zbin_offset=0.03, color1='red', color2='black', boxwidth=0.05, **kwds):
    """Plot boxplots of color distribution as a function of redshift bin for data 
       in two different dataframes.
       
       @param ax               axis object to plot on
       @param df1              dataframe 1
       @param df2              dataframe 2
       @param color_columns    name of column in df1,df2 containing color (list of 2 items) 
       @param redshift_columns name of column in df1,df2 containing redshifts (list of 2 items) 
       @param zmin             minimum redshift to start bins at
       @param zmax             maximum redshift to end bins at
       @param dbin             redshift bin width
       @param whis             see http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.boxplot
       @param zbin_offset      offset in redshift to plot second set of boxes at
       @param color1           color to plot boxplots from dataframe 1 (will be filled in)
       @param color2           color to plot boxplots from dataframe 2 (will just be outline)
       @param boxwidth         width (in redshift) to make boxplots
       @param kwds             other plotting keyword args to be passed to matplotlib boxplot func

    """
    
    # do redshift binning to start
    z_column1 = redshift_columns[0]
    z_column2 = redshift_columns[1]
    df1["cut_output"] = pd.cut(df1[z_column1], np.arange(zmin, zmax+dbin, dbin))
    df2["cut_output"] = pd.cut(df2[z_column2], np.arange(zmin, zmax+dbin, dbin))
    
    # redshift bin centers
    zbins = np.arange(zmin, zmax, dbin) + dbin/2.
    
    # color columns
    color_column1 = color_columns[0]
    color_column2 = color_columns[1]
    
    # plot the boxes from the first dataframe
    bp = df1.boxplot(column=[color_column1], by="cut_output", return_type='dict', 
                     ax=ax, positions=zbins, whis=whis, widths=boxwidth, **kwds)
                    
    # these boxes will be filled with 'color1'
    numBoxes = len(zbins)
    fill_boxes_with_color(numBoxes, bp, ax, color_column1, color1)

    
    # plot the boxes from the second dataframe
    bp = df2.boxplot(column=[color_column2], by="cut_output", return_type='dict', 
                    ax=ax, positions=zbins-zbin_offset, whis=whis, widths=boxwidth, **kwds)
              
    # these boxes will be outlined by 'color2'
    color_outlines(bp, color_column2, color2)
                    
    # set the tick labels along the redshift axes
    xtickNames = plt.setp(ax, xticklabels=zbins)
    
    # handmake the legend
    import matplotlib.lines as mpll

    xdata = np.arange(0,3,0.1)
    ydata = np.arange(0,3,0.1) 
    
    x1 = mpll.Line2D(xdata, ydata, linewidth=None, linestyle='none', color=color1, 
                marker='s', markersize=15, markeredgewidth=None, 
                markeredgecolor=None, markerfacecolor=None, 
                markerfacecoloralt='none', fillstyle=None, antialiased=None, 
                dash_capstyle=None, solid_capstyle=None, dash_joinstyle=None, 
                solid_joinstyle=None, pickradius=5, drawstyle=None, markevery=10,
                alpha=0.5)

    x2 = mpll.Line2D(xdata, ydata, linewidth=None, linestyle='solid', color=color2, 
                marker=None, markersize=15, markeredgewidth=None, 
                markeredgecolor=None, markerfacecolor=None, 
                markerfacecoloralt='none', fillstyle=None, antialiased=None, 
                dash_capstyle=None, solid_capstyle=None, dash_joinstyle=None, 
                solid_joinstyle=None, pickradius=5, drawstyle=None, markevery=10,
                alpha=1)
    ax.legend([x1, x2], ['data1', 'data2'])
    
    return [x1, x2]
