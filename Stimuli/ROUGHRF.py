from psychopy import visual, core, event, misc, logging, gui, monitors
import csv, pylab, math, random, numpy, time, imp, sys
sys.path.append("C:/Users/fitzlab1/Desktop/psychopy/triggers") #path to trigger classes
import serialTriggerDaqOut 

#---------- Stimulus Description -----------#
#info= {'stimType':1, 'logPath':'c:/users/fitzlab1/documents/', 'intTest':'exp001.txt'}
#infoDlg = gui.DlgFromDict(dictionary=info, title='Experiment Para
#---------- Monitor Properties ----------#
mon= monitors.Monitor('dualmonitors') #gets the calibration for stimMonitor 
#----------distance off by 2x, 50 cm == 25 cm
info={'stimType' :1}

#----------- Window Config -----------#
myWin = visual.Window(monitor='dualmonitors', 
    size=(3840,1080),
    fullscr=True,
    allowGUI = False,
    screen=1,
    units = 'deg') 

#INITIALISE SOME STIMULI
temporalFreq = 40 # temporal frequency in cycles/sec
#phaseShiftPeriod = 1# phase shift period for checkerboard
#phaseShiftPeriod1 = 2# phase shift period for stim
increment = 2.25	# initial increment for moving szddstims

myMask = numpy.array([
        [ -1,-1,1,-1],
        [ -1,-1,1,-1],
        [-1,-1,1, -1],
        [ -1,-1,1,-1],
        ])

#-----------Begin STIMULI -----------#
#----------- Gabor -----------#
if info['stimType'] ==1:
    stim = visual.PatchStim(myWin,pos=(-68,0), units = 'deg',
                           tex="sqr",mask="circle",
                           size=(60,60), sf=(0.05), ori=22.5*0,  depth=0, opacity=1.0,
                           autoLog=False)#this stim changes too much for autologging to work 'sqr', 'saw', 'tri'
#----------- End STIMULI -----------#           

stim2 = visual.PatchStim(myWin,pos=(68,-10), units = 'deg',
                           tex="sqr",mask=myMask,
                           size=(20,20), sf=(0.1), ori=22.5*0, depth=0, opacity=0,
                           autoLog=False)#this stim changes too much for autologging to work

#----------- Mouse setup -----------#ssssssswwwd
myMouse = event.Mouse(win=myWin)

count=0 # not sure what this is used for
contrastCounter=1;
#----------- Keyboard Setup -----------#
while True: #continue until keypress
    #handle key presses each frame
    for key in event.getKeys():
        #----------- All Stimuli -----------#
        if key in ['escape']:
            print "Stim Type", str(info['stimType'])
            print "Position", str(stim.pos)
            print "Orientation", str(stim.ori)
            print "SF", str(stim.sf)
            print "Size", str(stim.size)
            core.quit()
        if info['stimType'] ==1:
            if key in ['q']:
                stim.setSF(0.02, '+')
                print "SF = ", str(stim.sf)
                stim2.setSF(0.02, '+')
                print "SF = ", str(stim2.sf)
            if key in ['a']:
                stim.setSF(0.02, '-')
                print "SF = ", str(stim.sf)
                stim2.setSF(0.02, '-')
                print "SF = ", str(stim2.sf)
            if key in ['1']:
                stim.setSF(.08)
                print "SF = ", str(stim.sf)
                stim2.setSF(.08)
                print "SF = ", str(stim2.sf)
            if key in ['z']:
                contrastCounter=contrastCounter-1
                print "Stim is off"
                stim.setContrast(contrastCounter)
                stim2.setContrast(contrastCounter)
                if contrastCounter < 0:
                    contrastCounter=1
                    print "Stim is on"
                    stim.setContrast(contrastCounter)
                    stim2.setContrast(contrastCounter)
            if key in ['up']:
                stim.setPos([0,increment],'+')
                print "Position = ", str(stim.pos)
            if key in ['down']:
                stim.setPos([0,increment],'-')
                print "Position = ", str(stim.pos)
            if key in ['left']:
                stim.setPos([increment,0],'-')
                print "Position = ", str(stim.pos)
            if key in ['right']:
                stim.setPos([increment,0],'+')
                print "Position = ", str(stim.pos)
            if key in['8']:
                stim.setPos([0,0])
                print "Position=", str(stim.pos)
            if key in ['r']:
                temporalFreq += 1
                print "TF=", str(temporalFreq)
            if key in ['f']:
                temporalFreq -= 1
                print "TF=", str(temporalFreq)
            if key in ['4']:
                temporalFreq = 4
                print "TF=", str(temporalFreq)
            if key in ['t']:
                increment = increment + 0.5
                print "Position increment=", str(increment)
            if key in ['g']:
                increment = increment  - 0.5
                print "Position increment=", str(increment)
            if key in ['5']:
                increment = 0.5
                print "Position increment=", str(increment)
            if key in ['e']:
                stim.setOri(11.25, '+')
                print "Orientation=", str(stim.ori+90)
            if key in ['d']:
                stim.setOri(11.25, '-')
                print "Orientation=", str(stim.ori+90)
            if key in['c']:
                stim.setOri(90,'-')
                print "Orientation=", str(stim.ori+90)
            if key in ['3']:
                stim.setOri(-90)
                print "Orientation=", str(stim.ori+90)
            if key in ['s']:
                stim.setSize(1,'-')
                print "Size = ", str(stim.size)
                stim2.setSize(1,'-')
                print "Size = ", str(stim2.size)
            if key in ['x']:
                stim.setSize(500)
                print "Size = ", str(stim.size)
            if key in ['2']:
                stim.setSize(10)
                print "Size = ", str(stim.size)
                stim2.setSize(10)
                print "Size = ", str(stim2.size)
            if key in ['w']:
                stim.setSize(1, '+')
                print "Size = ", str(stim.size)
                stim2.setSize(1, '+')
                print "Size = ", str(stim2.size)
                
            if key in ['b']:
                stim.setPhase(0.05, '+')
                print "Phase=", str(stim.phase)
            if key in ['n']:
                stim.setPhase(0.05, '-')
                print "Phase=", str(stim.phase)
            if key in ['m']:
                stim.setPhase(0)
                print "Phase=", str(stim.phase)
                
            if key in ['y']:
                stim2.setOri(11.25, '+')
                print "Orientation=", str(stim2.ori+90)
            if key in ['h']:
                stim2.setOri(11.25, '-')
                print "Orientation=", str(stim2.ori+90)
            if key in ['6']:
                stim2.setOri(-90)
                print "Orientation=", str(stim2.ori+90)
            if key in ['9']:
                stim2.setPos([0,increment],'+')
                print "Position = ", str(stim2.pos)
            if key in ['o']:
                stim2.setPos([0,increment],'-')
                print "Position = ", str(stim2.pos)
            if key in ['i']:
                stim2.setPos([increment,0],'-')
                print "Position = ", str(stim2.pos)
            if key in ['p']:
                stim2.setPos([increment,0],'+')
                print "Position = ", str(stim2.pos)
            if key in['0']:
                stim2.setPos([0,0])
                print "Position=", str(stim2.pos)
            if key in ['j']:
                stim2.setPhase(0.05, '+')
                print "Phase=", str(stim2.phase)
            if key in ['k']:
                stim2.setPhase(0.05, '-')
                print "Phase=", str(stim2.phase)
            if key in ['l']:
                stim2.setPhase(0)
                print "Phase=", str(stim2.phase)
        #----------- Checkerboard Stimulus -----------#

    #get mouse events
    #Handle the wheel(s):
    # Y is the normal mouse wheel, but some (e.g. mighty mouse) have an x as well
    mouse_dX,mouse_dY = myMouse.getRel()
    mouse1, mouse2, mouse3 = myMouse.getPressed()
    if (mouse3):
        stim.setSize(mouse_dX, '+',1)
        stim.setPos([mouse_dX, mouse_dY], '+', 20)
    wheel_dX, wheel_dY = myMouse.getWheelRel()
    stim.setOri(wheel_dY*5, '+')
    clock = core.Clock()
    event.clearEvents()#get rid of other, unprocessed events
    #do the draw
    if info['stimType'] == 1:
        stim.setPhase(.001*temporalFreq,'+')
        stim.draw()
        myWin.flip()#redraw the buffer
        
#        stim2.setPhase(.001*temporalFreq,'+')
#        stim2.draw()
#        myWin.flip()#redraw the buffer

    else:
        stim.draw()
        myWin.flip()
        count += 1
