from psychopy import visual, logging, core, monitors
import math, random, numpy, sys

sys.path.append("Z:/Juliane/Code/Stimuli/triggers")  # path to trigger classes
from os import path

print "initialized"

# ---------- Stimulus Description ---------- #
'''A fullscreen drifting grating for 2pt orientation tuning'''
#Experiment logging parameters
dataPath='x:/'
animalName='F2267_2018-10-22'; 
logFilePath =dataPath+animalName+'\\'+animalName+'.txt' #including filepath

#---------- Monitor Properties ----------#
mon= monitors.Monitor('testMonitor') #gets the calibration for stimMonitor
mon.setDistance(30)
overwriteGammaCalibration = False
newGamma = 0.479
myWin = visual.Window(size=[1920, 1080], monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=False, units= 'pxl')
# ---------- Stimulus Parameters ---------- #
# trials and duration
numDirections = 8  # typically 4, 8, or 16#
# orientations = numpy.arange(0.0,180,180.0/numOrientations) #Remember, ranges in Python do NOT include the final value!
directions = numpy.arange(0, 360,
                            360.0 / numDirections)  # Remember, ranges in Python do NOT include the final value!
# orientations = numpy.arange(90,91,1);
# orientations = [90,270]

numTrials = 6  # Run all the stims this many times

doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.
nBlank = 1  # number of blanks to show per trial.

stimDuration = 3
isi =3

isRandom = 1
initialDelay = .5  # time in seconds to wait before first stimuli. Set to 0 to begin ASAP.

# Grating parameter
contrast = 1
# aperture and position parameters
# Triggering type
# Can be any of:#"NoTrigger" - no triggering; stim will run freely
# "SerialDaqOut" - Triggering by serial port. Stim codes are written to the MCC DAQ.
# "OutOnly" - no input trigger, but does all output (to CED) and logging
# "DaqIntrinsicTrigger" - waits for stimcodes on the MCC DAQ and displays the appropriate stim ID
# triggerType = 'DaqIntrinsicTrigger'
triggerType = 'None'
serialPortName = 'COM3'  # ignored if triggerType is "None"
adjustDurationToMatch2P = True
# Experiment logging parameters
import time

date = (time.strftime("%Y-%m-%d"))

# ---------- Stimulus code begins here ---------- #
stimCodeName = path.dirname(path.realpath(__file__)) + '\\' + path.basename(__file__)

# make a window
if overwriteGammaCalibration:
    myWin.setGamma(newGamma)
    print "Overwriting Gamma Calibration. New Gamma value:", newGamma

print "made window, setting up triggers"
#Set up the trigger behavior
trigger = 'NoTrigger'
if triggerType == "NoTrigger":
    import noTrigger
    trigger = noTrigger.noTrigger(None) 
elif triggerType == "SerialDaqOut" or triggerType == 'OutOnly':
    import serialTriggerDaqOut
    print 'Imported trigger serialTriggerDaqOut'
    trigger = serialTriggerDaqOut.serialTriggerDaqOut(serialPortName) 
    # determine the Next experiment file name
    expName=trigger.getNextExpName([dataPath,animalName])
    print "Trial name: ",expName
    if triggerType == 'OutOnly':
        trigger.readSer=False
    #Record a bunch of serial triggers and fit the stim duration to an exact multiple of the trigger time
    if adjustDurationToMatch2P:
        print "Waiting for serial Triggers"
    # store the stimulus data and prepare the directory
    trigger.preTrialLogging([dataPath,animalName,expName,stimCodeName,directions,logFilePath])
elif triggerType=="DaqIntrinsicTrigger":
    import daqIntrinsicTrigger
    trigger = daqIntrinsicTrigger.daqIntrinsicTrigger(None) 
else:
    print "Unknown trigger type", triggerType
        
print stimDuration

# create grating stim
dotPatch = visual.DotStim(win=myWin, color = (-1.0, -1.0, -1.0), units = 'deg', contrast = 1, dir = 0, nDots= 300,
                          dotSize= 15, fieldShape= 'circle', fieldSize= 50)
dotPatch.setAutoDraw(True)
dotPatch.mask = 'circle'

clrctr = 1
print "made grating"
# run
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
print "\n", str(len(directions) + doBlank), "stims will be run for", str(numTrials), "trials."
if nBlank > 1:
    print "Will run blank " + str(nBlank) + " times"
duration = initialDelay + numTrials * numDirections * (isi + stimDuration) + numTrials * nBlank
print"estimated duration: " + str(math.ceil(duration / 60)) + " min/" + str(duration) + " seconds"
# force a wait period of at least 5 seconds before first stim
if initialDelay > 0:
    print" waiting " + str(initialDelay) + " seconds before starting stim to acquire a baseline."
    dotPatch.setContrast(0)
    # flipStim.setContrast(0)
    while clock.getTime() < initialDelay:
        myWin.flip()
for trial in range(0, numTrials):
    # determine stim order
    print "Beginning Trial", trial + 1
    stimOrder = range(0, len(directions) + doBlank)
    if nBlank > 1:
        blankID = len(directions)
        for ibl in range(1, nBlank):
            stimOrder.append(blankID)
    if isRandom:
        random.shuffle(stimOrder)
    for stimNumber in stimOrder:
        # display each stim
        #trigger.preStim(stimNumber + 1)
        # display stim
        # flipStim.setContrast(1)
        # flipStim.setAutoDraw(True)
        # convert orientations to standard lab notation
        if stimNumber == len(directions):
            dotPatch.setContrast(0)
            print "\tStim", stimNumber + 1, " (blank)"
        else:
            dotPatch.setContrast(contrast)
            dotPatch.dir = directions[stimNumber] - 90
            print "\tStim", stimNumber + 1, directions[stimNumber], 'deg'
        clock.reset()
        while clock.getTime() < stimDuration:
            clrctr = clrctr + 1;
            # if clrctr%2==1:
            # flipStim.setColor((0,0,0),colorSpace='rgb')
            # flipStim.setContrast(-1)
            # else:
            # flipStim.setColor((1,1,1),colorSpace='rgb')
            # flipStim.setContrast(1)
            #trigger.preFlip(None)
            myWin.flip()
            #trigger.postFlip(None)

        # now do ISI
        if isi != 0:
            clock.reset()
            print 'ISI'
            dotPatch.setContrast(0)
            # flipStim.setContrast(0)
            #trigger.preFlip(None)
            myWin.flip()
            #trigger.postFlip(None)
            # flipStim.setAutoDraw(False)
            while clock.getTime() < isi:
                # Keep flipping during the ISI. If you don't do this you can get weird flash artifacts when you resume flipping later.
                myWin.flip()
        #trigger.postStim(None)

#trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'