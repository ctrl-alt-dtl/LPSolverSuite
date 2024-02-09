#/------------------------------------------------------------/#
#Name: Nash Arbitration Solver v1.2
#Creator: David "Ty" Long
#Created: 11 Dec 2016
#Last Edit: 17 Oct 2022
#Dependencies: Pyhton 2.7 (Updated for Python 3.11 10/17/2022),
#              Matplotlib, & Numpy
#            Use 'easy_install' or 'pip' to download install or
#            Download Matplotlib from: http://matplotlib.org/
#            Download Numpy from: http://www.numpy.org/
#            
#            pip3.11 install matplotlib, numpy, tk
#            Homebrew: brew install python-tk@3.11
#/------------------------------------------------------------/#

from __future__ import division
import tkinter as tk
from tkinter import StringVar 
import sys
from fractions import Fraction
import random
from tkinter import messagebox
import matplotlib
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
import numpy as np
from numpy import vstack, ones
from numpy.linalg import lstsq

#import matplotlib.pyplot as plt # DEAD CODE
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

CLOCKWISE = -1
COLLINEAR = 0
COUNTERCLOCKWISE = +1
eps = sys.float_info.epsilon

class HelpView(tk.Frame):
    #Calls using two, only need self here. 
    #@TODO: Fix
    
    #############################################
    #
    #############################################
    def showHelp(self, dumbArg): 
        tk.Frame.__init__(self)
        t = tk.Toplevel(self)
        t.wm_title("Help Menu")
        l = tk.Label(t, text="Help String here.")
        l.pack(side="top",fill="both", expand=True)

class MainView(tk.Frame):
    #Globals
    global cellA, cellB, cellC, cellD, sqCell1, sqCell2
    global wMethodAB, wMethodBA, wMethodCD, wMethodDC
    global deltaA, deltaB, deltaC, deltaD
    global cellA1, cellA2, cellB1, cellB2
    global cellC1, cellC2, cellD1, cellD2
    global rTotal, cTotal
    
    #Initialization
    wMethodAB = None
    wMethodBA = None
    wMethodCD = None
    wMethodDC = None
    deltaA = deltaB = deltaC = deltaD = None
    cellA1 = cellA2 = cellB1 = cellB2 = None
    cellC1 = cellC2 = cellD1 = cellD2 = None
    sqCell1 = sqCell2 = None
    rTotal = cTotal = None
    
    def __init__(self, *args, **kwargs):
        '''
        Non-tested code
        # valid percent substitutions (from the Tk entry man page)
        # note: you only have to register the ones you need; this
        # example registers them all for illustrative purposes
        #
        # %d = Type of action (1=insert, 0=delete, -1 for others)
        # %i = index of char string to be inserted/deleted, or -1
        # %P = value of the entry if the edit is allowed
        # %s = value of entry prior to editing
        # %S = the text string being inserted or deleted, if any
        # %v = the type of validation that is currently set
        # %V = the type of validation that triggered the callback
        #      (key, focusin, focusout, forced)
        # %W = the tk name of the widget
        
        def onValidate(self, action, index, value_if_allowed,
                       prior_value, text, validation_type, trigger_type, widget_name):
            if text in '0123456789.':
                try:
                    float(value_if_allowed)
                    return True
                except ValueError:
                    return False
            else:
                return False
        '''
            
        #StringVar Initializations
        wMethodAB = StringVar()
        wMethodBA = StringVar()
        wMethodCD = StringVar()
        wMethodDC = StringVar()
        deltaA = StringVar()
        deltaB = StringVar()
        deltaC = StringVar()
        deltaD = StringVar()
        cellA1 = StringVar()
        cellA2 = StringVar()
        cellB1 = StringVar()
        cellB2 = StringVar()
        cellC1 = StringVar()
        cellC2 = StringVar()
        cellD1 = StringVar()
        cellD2 = StringVar()
        rTotal = StringVar()
        cTotal = StringVar()
        sqCell1 = StringVar()
        sqCell2 = StringVar()
        
        #############################################
        #Function used for Gift Wrapping (a.k.a - Jarvis March) Algorithm
        #############################################
        def orientation(a, b):
            x0, y0 = a
            x1, y1 = b
            cross = x0 * y1 - x1 * y0
            if cross > eps:
                return COUNTERCLOCKWISE
            elif cross < -eps:
                return CLOCKWISE
            else:
                return COLLINEAR
        
        #############################################
        #Function used for Gift Wrapping (a.k.a - Jarvis March) Algorithm
        #############################################
        def same_halfplane(a, b):
            x0, y0 = a
            x1, y1 = b
            dot = x0 * x1 + y0 * y1
            if dot >= eps:
                return True
            elif dot < eps:
                return False
        
        #############################################
        #Gift Wrapping (a.k.a - Jarvis March) Algorithm
        #############################################
        def jarvis(points):
            points = points[:]
            r0 = min(points)
            hull = [r0]
            r, u = r0, None
            remainingPoints = [x for x in points if x not in hull]
            while u != r0 and remainingPoints:
                u = random.choice(remainingPoints)
                for t in points:
                    a = (u[0] - r[0], u[1] - r[1])
                    b = (t[0] - u[0], t[1] - u[1])
                    if (t != u and
                        (orientation(a, b) == CLOCKWISE or
                         (orientation(a, b) == COLLINEAR and
                          same_halfplane(a, b)))):
                        u = t
                r = u
                points.remove(r)
                hull.append(r)
                try:
                    remainingPoints.remove(r)
                except ValueError:
                    # ValueError: list.remove(x): x not in list
                    pass
            return hull
                
        #############################################
        #Locate the distance of points to determine line intersections
        #############################################
        def linePoints(self, pSet1, pSet2):
            A = (pSet1[1] - pSet2[1])
            B = (pSet2[0] - pSet1[0])
            C = (pSet1[0]*pSet2[1] - pSet2[0]*pSet1[1])
            return A, B, -C
        
        #############################################
        #Calculate the intersections of two lines
        #############################################
        def lineIntersection(self, lSet1, lSet2):
            D  = lSet1[0] * lSet2[1] - lSet1[1] * lSet2[0]
            Dx = lSet1[2] * lSet2[1] - lSet1[1] * lSet2[2]
            Dy = lSet1[0] * lSet2[2] - lSet1[2] * lSet2[0]
            if D != 0:
                x = Dx / D
                y = Dy / D
                return x,y
            else:
                return False
        
        #############################################
        #Find the equation of a line y = mx + b
        #############################################
        def lineEquation(self, x_coords, y_coords):
            #x_coords, y_coords = zip(*points)
            A = vstack([x_coords, ones(len(x_coords))]).T
            m, b = lstsq(A, y_coords)[0]
            m=str(Fraction(m).limit_denominator())
            #DEBUG
            print("The equation of a line is:\ny = {m}x + {b:.2f}".format(m=m, b=b))
            lineEq = "y = {m}x + {b:.2f}".format(m=m, b=b)
            tk.Label(self, text=str(lineEq)).grid(row=6, column=12, ipadx=7)
            
        #############################################
        #Displays all plot points and lines on the matplotlib graph
        #############################################
        def plotPoints_Lines(self, a1, a2, b1, b2, \
                             c1, c2, d1, d2, sqRow, sqCol):
            #Globals
            nfpCoord1 = None
            nfpCoord2 = None
            
            #Plot setup
            maxValX = max(a1, a2, b1, b2)
            maxValY = max(c1, c2, d1, d2)
            minValX = min(a1, a2, b1, b2)
            minValY = min(c1, c2, d1, d2)
            
            #Plot area extensions
            maxValX += 5
            maxValY += 5
            minValX -= 5
            minValY -= 5
            
            #Creating SQ line plot segment
            sqLineX = [(sqRow, sqCol), (maxValX, sqCol)]
            sqLineY = [(sqRow, sqCol), (sqRow, maxValY)]
            (line1_xs, line1_ys) = list(zip(*sqLineX))
            (line2_xs, line2_ys) = list(zip(*sqLineY))
            
            #Embedded Plot Graph
            points = ([a1, a2, b1, b2], [c1, d1, c2, d2])
            points = list(zip(*points))
            hull = jarvis(points)
            hx, hy = list(zip(*hull))
            #The points of the last two locations counterclockwise
            fhx, fhy = zip(hull[3], hull[0])
            
            #debug Convex Hull
            print("\nHull Order:")
            print(list(hull))
            
            #Sort the points to find the highest two
            max_x = (sorted(points, key=lambda x: x[0]))[-1]
            max_y = (sorted(points, key=lambda x: x[1]))[-1]
            
            #Find the equation of a line between two points
            lineEquation(self, max_x, max_y)
            
            #Extend the Pareto Optimal Line to the SQ Level Set
            pSet1 = ([sqRow, sqCol]) 
            pSet2 = ([sqRow, maxValY])
            lineXset = linePoints(self, pSet1, pSet2)
            pSet1 = ([sqRow, sqCol]) 
            pSet2 = ([maxValX, sqCol])
            lineYset = linePoints(self, pSet1, pSet2)
            paretoOptimalPoints = linePoints(self, max_x, max_y)
            intersectX = lineIntersection(self, paretoOptimalPoints, lineXset)
            intersectY = lineIntersection(self, paretoOptimalPoints, lineYset)
            try:
                po_x1, po_y1 = list(zip(intersectX))
                po_x2, po_y2 = list(zip(intersectY))
            except TypeError:
                messagebox.showerror("Error", \
                                               "Values Incorrect. " + \
                                               "Check Inputs!")
                print("Invalid Inputs. Iteration not possible.")
            pOptimalSet1, pOptimalSet2 = list(zip(intersectX,intersectY))
            
            #Map the largest point on both Row X and Column Y axis
            #To use as a comparison tool for determining NFP
            x1, y1= map(list, list(zip(max_x)))
            x2, y2 = map(list, list(zip(max_y)))
            
            #Negotiation Set
            if(sqRow > np.array([x1])):
                nfpCoord1 = max_x
            else:
                nfpCoord1 = intersectX
                
            if(sqCol > np.array([y2])):
                nfpCoord2 = max_y
            else:
                nfpCoord2 = intersectY
            
            nashFairPoint = nfpSolver(self, nfpCoord1, nfpCoord2)
            nfp_x, nfp_y = list(zip(*nashFairPoint))
            
            #Configure Widget
            fig = Figure(figsize=(5,5), dpi=100)
            ax = fig.add_subplot(1,1,1)
            ax.grid(which='both')
            #Plot Points
            ax.plot(a1, c1, 'ko')
            ax.plot(a2, d1, 'ko')
            ax.plot(b1, c2, 'ko')
            ax.plot(b2, d2, 'ko')
            #SQ Plot
            ax.plot(sqRow, sqCol, 'ro')
            
            #Extend Pareto Optimal lines to SQ lines
            ax.plot(pOptimalSet1, pOptimalSet2, 'k--')
            
            #Convex Hull Plot
            ax.plot(hx, hy, 'k.-', markersize=10)
            #Close polygon
            ax.plot(fhx, fhy, 'k.-', markersize=10)
            #Label Points
            ax.text(a1+.10, c1+.10, "(A,C)", fontsize=8)
            ax.text(a1+.10, c1+.10, "(A,C)", fontsize=8)
            ax.text(a2+.10, d1+.10, "(A,D)", fontsize=8)
            ax.text(b1+.10, c2+.10, "(B,C)", fontsize=8)
            ax.text(b2+.10, d2+.10, "(B,D)", fontsize=8)
            ax.text(sqRow-1, sqCol, "(SQ)", fontsize=8)

            #Draw SQ Lines
            ax.add_line(Line2D(line1_xs, line1_ys, linewidth=2, color='green'))
            ax.add_line(Line2D(line2_xs, line2_ys, linewidth=2, color='red'))
            
            #Pareto Optimal Plot
            ax.plot(po_x1, po_y1, 'bo')
            ax.plot(po_x2, po_y2, 'bo')
            
            #Nash Fair Point
            ax.plot(nfp_x, nfp_y, 'go')
            ax.text(np.add(nfp_x,.2), np.add(nfp_y,.2), "NFP", fontsize=8)
            nfpText = "NFP: {}".format(str(nashFairPoint).strip('[]'))
            tk.Label(self, text=str(nfpText)).grid(row=7, column=12, ipadx=12)
            
            #Extend Plot Area 
            ax.set_xbound(minValX, maxValX)
            ax.set_ybound(minValY, maxValY)
            #Draw Canvas
            canvas = FigureCanvasTkAgg(fig, self)
            canvas.get_tk_widget().grid(row=9, columnspan=15, column=0, sticky="s")
            canvas.show()
            
            #Crash!! (if resized)
            toolbar = NavigationToolbar2Tk(canvas, self)
            toolbar.update()
            toolbar.grid(row=10, columnspan=15, column=1, sticky="w")
            #Major Crash!!
            #canvas._tkcanvas.grid(row=10, columnspan=15, column=1, sticky="w")
            
            '''
            #NOT NEEDED (Dead Code)
            #PyPlot Pop-up Window
            fig2 = plt.figure()
            plt.grid(True)
            plt.title("Graph Plot")
            plt.xlabel("Row")
            plt.ylabel("Column")
            plt.plot([a1, a2, b1, b2], [c1, d1, c2, d2], 'ko')
            plt.plot([sqRow], [sqCol], 'ro')
            plt.plot(hx, hy, 'k.-', markersize=10)
            plt.plot(fhx, fhy, 'k.-', markersize=10)
            #Label points
            plt.text(a1+.10, c1+.10, "(A,C)")
            plt.text(a2+.10, d1+.10, "(A,D)")
            plt.text(b1+.10, c2+.10, "(B,C)")
            plt.text(b2+.10, d2+.10, "(B,D)")
            plt.text(sqRow-1, sqCol, "(SQ)")
            
            plt.axis([minValX, maxValX, minValY, maxValY], aspect=1)
            
            #@TODO: Fix plot lines for pyplot for when points are not 1:1
            #PyPlot for SQ Lines
            plt.axvline(sqRow, ymin=.5, ymax=1, color='red')
            plt.axhline(sqCol, xmin=.5, xmax=1, color='green')
            
            fig2.canvas.set_window_title("Nash Arbitration Plots")
            plt.show()
            '''
        
        #############################################
        #Solve for Nash Fair Point along the negotiation 
        #############################################
        def nfpSolver(self, nfpCoord1, nfpCoord2):
            x1, y1 = map(list, zip(nfpCoord1))
            x2, y2 = map(list, zip(nfpCoord2))
            x_m_point = np.add(x1,x2)/2
            y_m_point = np.add(y1,y2)/2
            round_x_point = [round(xprec, 2) for xprec in x_m_point]
            round_y_point = [round(yprec, 2) for yprec in y_m_point]
            nashFairPoint = list(zip(round_x_point, round_y_point))
            return nashFairPoint
        
        #############################################
        #Solves for the Security Level
        #############################################
        def sqSolver(self, sqOverRow, sqOverCol):
            #Function Variables
            sqRow = None
            sqCol = None
            #subscript_0 = u'\u2080' #Future implementation of subscript 0
            a1 = int(cellA1.get())
            a2 = int(cellA2.get())
            b1 = int(cellB1.get())
            b2 = int(cellB2.get())
            c1 = int(cellC1.get())
            c2 = int(cellC2.get())
            d1 = int(cellD1.get())
            d2 = int(cellD2.get())
            
            if not sqOverRow:
                #Security Level Row Side
                #Row is minimizing and Column is maximizing
                sqRowAmin = min(a1, a2)
                sqRowBmin = min(b1, b2)
                sqColCmax = max(a1, b1)
                sqColDmax = max(a2, b2)
                
                if(sqRowAmin == sqColCmax):
                    sqRow = sqRowAmin
                    rowSqText = "X0: {:d}".format(sqRow)
                elif(sqRowAmin == sqColDmax):
                    sqRow = sqRowAmin
                    rowSqText = "X0: {:d}".format(sqRow)
                elif(sqRowBmin == sqColCmax):
                    sqRow = sqRowBmin
                    rowSqText = "X0: {:d}".format(sqRow)
                elif(sqRowBmin == sqColDmax):
                    sqRow = sqRowBmin
                    rowSqText = "X0: {:d}".format(sqRow)
                #No Prudential Strategy, move on to Expected Values
                else:
                    #Expected Values Row Side
                    sqRowEv1 = a1*(float(Fraction(wMethodAB.get()))) \
                            + b1*(float(Fraction(wMethodBA.get())))
                    sqRowEv2 = a2*(float(Fraction(wMethodAB.get()))) \
                            + b2*(float(Fraction(wMethodBA.get())))
                    
                    #Mainly for precision comparison since similar float variables \
                    #can be different at higher precisions. We only need 0.000     \
                    #in precision, but this doesn't hurt.                          \
                    rowAbs = abs(sqRowEv1 - sqRowEv2)
                    if rowAbs < 0.000001:
                        sqRow = sqRowEv1
                        rowSqText = "X0: {:.2f}".format(sqRow)
                    else:
                        #Unlikely this will ever pop up
                        messagebox.showerror("Error", \
                                               "Expected Values Do Not Match.\n" + \
                                               " ABS(EV(A+B) is > %d" % rowAbs )
            else:
                rowSqText = "X0: {:.2f}".format(sqOverRow)
                sqRow = sqOverRow
            
            if not sqOverCol:
                #Security Level Row Side
                #Row is maximizing and Column is minimizing
                sqColCmin = min(c1, c2)
                sqColDmin = min(d1, d2)
                sqRowAmax = max(c1, d1)
                sqRowBmax = max(c2, d2)
                
                if(sqColCmin == sqRowAmax):
                    sqCol = sqColCmin
                    colSqText = "Y0: {:d}".format(sqCol)
                elif(sqColCmin == sqRowBmax):
                    sqCol = sqColCmin
                    colSqText = "Y0: {:d}".format(sqCol)
                elif(sqColDmin == sqRowAmax):
                    sqCol = sqColDmin
                    colSqText = "Y0: {:d}".format(sqCol)
                elif(sqColDmin == sqRowBmax):
                    sqCol = sqColDmin
                    colSqText = "Y0: {:d}".format(sqCol)
                #No Prudential Strategy, move on to Expected Values
                else:
                    #Expected Values Column Side
                    sqColEv1 = c1*(float(Fraction(wMethodCD.get()))) \
                            + d1*(float(Fraction(wMethodDC.get())))
                    sqColEv2 = c2*(float(Fraction(wMethodCD.get()))) \
                            + d2*(float(Fraction(wMethodDC.get())))
                    
                    #Mainly for precision comparison since similar float variables \
                    #can be different at higher precisions. We only need 0.000     \
                    #in precision, but this doesn't hurt.                          \
                    colAbs = abs(sqColEv1 - sqColEv2)
                    if colAbs < 0.000001:
                        sqCol = sqColEv1
                        colSqText = "Y0: {:.2f}".format(sqCol)
                    else:
                        #Unlikely this will ever pop up
                        messagebox.showerror("Error", \
                                               "Expected Values Do Not Match.\n" + \
                                               " ABS(EV(C+D) is > %d" % colAbs )
            else:
                colSqText = "Y0: {:.2f}".format(sqOverCol)
                sqCol = sqOverCol
                    
            tk.Label(self, text=str(rowSqText)).grid(row=4, column=12, ipadx=7)
            tk.Label(self, text=str(colSqText)).grid(row=5, column=12, ipadx=7)
            
            #Send the calculated data to plot
            plotPoints_Lines(self, a1, a2, b1, b2, c1, \
                             c2, d1, d2, sqRow, sqCol)
        
        #############################################
        #Overriding Security Level
        #############################################
        def sqOverride(self, sqOverRow, sqOverCol):
            sqOverRow = int(sqCell1.get())
            sqOverCol = int(sqCell2.get())
            sqSolver(self, sqOverRow, sqOverCol)
        
        #############################################
        #Solves for Williams Method of play Strategy
        #############################################
        def williamsSolver(self, cellA1, cellA2, cellB1, cellB2, \
                           cellC1, cellC2, cellD1, cellD2):
            
            deltaA.set(abs(abs(int(cellA1.get())) - abs(int(cellA2.get()))))
            deltaB.set(abs(abs(int(cellB1.get())) - abs(int(cellB2.get()))))
            deltaC.set(abs(abs(int(cellC1.get())) - abs(int(cellC2.get()))))
            deltaD.set(abs(abs(int(cellD1.get())) - abs(int(cellD2.get()))))
            
            rTotal.set(str(int(deltaA.get())+int(deltaB.get())))
            cTotal.set(str(int(deltaC.get())+int(deltaD.get())))
            
            wMethodAB.set(str(Fraction(int(deltaB.get()), int(rTotal.get()))))
            wMethodBA.set(str(Fraction(int(deltaA.get()), int(rTotal.get()))))
            wMethodCD.set(str(Fraction(int(deltaD.get()), int(cTotal.get()))))
            wMethodDC.set(str(Fraction(int(deltaC.get()), int(cTotal.get()))))
        
        
        ########
        # Init Display Frame
        ########
        
        #Frame Initialization
        tk.Frame.__init__(self, *args, **kwargs)
        
        
        #Form setup and stucture
        tk.Label(self, text="Column", bg="red").grid(row=0, column=4)
        tk.Label(self, text="Row", bg="green").grid(row=3, column=0)
        tk.Label(self, text="C").grid(row=1, column=3)
        tk.Label(self, text="D").grid(row=1, column=6)
        tk.Label(self, text="A").grid(row=2, column=1)
        tk.Label(self, text="B").grid(row=4, column=1)
        
        #Additional labels
        tk.Label(self, text=u'\u0394'+"  ", anchor="w").grid(row=1, column=8)
        tk.Label(self, text=u'\u0394').grid(row=5, column=1)
        tk.Label(self, text="T=").grid(row=6, column=1)
        tk.Label(self, text="T=").grid(row=6, column=7)
        
        #Use for Delta change calculation
        tk.Label(self, textvariable=deltaA).grid(row=2, column=8, sticky="w")
        tk.Label(self, textvariable=deltaB).grid(row=4, column=8, sticky="w")
        tk.Label(self, textvariable=deltaC).grid(row=5, column=3, sticky="w")
        tk.Label(self, textvariable=deltaD).grid(row=5, column=6, sticky="w")
        
        tk.Label(self, textvariable=rTotal, width=2, anchor="w").grid(row=6, column=8)
        tk.Label(self, textvariable=cTotal, anchor="w").grid(row=6, column=2)
        
        '''
        Non-Tested Code
        vcmd = (self.register(onValidate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        #Entry AC
        tk.Entry(self, width=3, textvariable=cellA1, validate='key', validatecommand=vcmd).grid(row=2, column=2, sticky="w")
        '''
        
        #Entry AC
        tk.Entry(self, width=3, textvariable=cellA1).grid(row=2, column=2, sticky="w")
        tk.Entry(self, width=3, textvariable=cellC1).grid(row=2, column=3, sticky="w")
        #Entry AD
        tk.Entry(self, width=3, textvariable=cellA2).grid(row=2, column=5, sticky="w")
        tk.Entry(self, width=3, textvariable=cellD1).grid(row=2, column=6, sticky="w")
        #Entry BC
        tk.Entry(self, width=3, textvariable=cellB1).grid(row=4, column=2, sticky="w")
        tk.Entry(self, width=3, textvariable=cellC2).grid(row=4, column=3, sticky="w")
        #Entry BD
        tk.Entry(self, width=3, textvariable=cellB2).grid(row=4, column=5, sticky="w")
        tk.Entry(self, width=3, textvariable=cellD2).grid(row=4, column=6, sticky="w")
        
        
        #Williams Method spaces
        tk.Label(self, width=4, bg="black", fg="white", textvariable=wMethodAB)\
        .grid(row=2, column=9, sticky="w")
        tk.Label(self, width=4, bg="black", fg="white", textvariable=wMethodBA)\
        .grid(row=4, column=9, sticky="w")
        tk.Label(self, width=4, bg="black", fg="white", textvariable=wMethodCD)\
        .grid(row=6, column=3, sticky="w")
        tk.Label(self, width=4, bg="black", fg="white", textvariable=wMethodDC)\
        .grid(row=6, column=6, sticky="w")
        
        #It wouldn't be a GUI without buttons
        tk.Button(self, text="Find Williams", fg="red", command=lambda: \
                  williamsSolver(self, cellA1, cellA2, cellB1, cellB2, \
                                 cellC1, cellC2, cellD1, cellD2))\
        .grid(row=1, ipadx=7, column=12)
        tk.Button(self, text="Solve", command=lambda: sqSolver(self, None, None)) \
        .grid(row=2, ipadx=7, column=12)
        tk.Entry(self, width=2, textvariable=sqCell1).grid(row=3, column=13, sticky="w")
        tk.Entry(self, width=2, textvariable=sqCell2).grid(row=3, column=14, sticky="w")
        tk.Button(self, text="SQ Override", command=lambda: sqOverride(self, sqCell1, sqCell2)) \
        .grid(row=3, column=12)
        
        
        
        #Canvas widget for displaying plots and data
        tk.Label(self, text="Graph Plot").grid(row=8, column=6)
        tk.Canvas(self, width=520, height=500, bg="gray").grid(row=9, columnspan=15, column=0, sticky="s")
        
if __name__ == "__main__":
    print("Initializing chache. Please do not close this window...\n")
    #Tkinter Initialization
    root = tk.Tk()
    
    #Menu Bar Setup
    menuBar = tk.Menu(root)
    fileMenu = tk.Menu(menuBar, tearoff=0)
    menuBar.add_cascade(label="File", menu=fileMenu)
    helpMenu = tk.Menu(menuBar, tearoff=0)
    menuBar.add_cascade(label="Help", menu=helpMenu)
    aboutMenu = tk.Menu(menuBar, tearoff=0)
    menuBar.add_cascade(label="About", menu=aboutMenu)
    root.config(menu=menuBar)
    
    #Main Window
    view = MainView(root)
    view.pack(side="top", fill="both", expand=True)
    root.wm_geometry("600x800")
    root.maxsize(1100, 900)
    root.title("Nash Arbitration Solver")
    root.bind("<F1>", HelpView().showHelp)
    root.mainloop()
