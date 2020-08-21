import serial, csv, time, math, datetime
from psychopy import visual, core, filters, monitors
import pylab, math, random, numpy, sys

sys.path.append("Z:/Juliane/Code/Stimuli/triggers")
from os import path

print "initialized"

# ---------- Stimulus Description ---------- #
triggerType = 'NoTrigger'  # 'SerialDaqOut' OR 'OutOnly'
serialPortName = 'COM2'  # ignored if triggerType is "COM3"
adjustDurationToMatch2P = True

# Experiment logging parameters
import time

dataPath = '//mpfiwdtfitz68/Spike2Data/'
date = (time.strftime("%Y-%m-%d"))
logFilePath = dataPath + date + '\\' + date + '.txt'  # including filepath
stimCodeName = path.dirname(path.realpath(__file__)) + '\\' + path.basename(__file__)

mon = monitors.Monitor('testMonitor')
mon.setDistance(25)

# make a window
myWin = visual.Window(size=[1920, 1080], monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=False)
print "made window"

#experiment parameters
doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.
changeDirectionAt = 1  # When do we change movement directions? If 1, there should be no reversal. Setting to 0 causes instantaneous reversal, effectively inverting all stimIDs and angles.
stimDuration = .25
isi = 1
isRandom = 1
initialDelay = 5
numTrials = 10

#grating parameters
numWedges = 6
stimSize = [60,60]
centerPoint = [0, 0]
temporalFreq = 3
spatialFreq = 0.1
contrast = 1
textureType = 'sqr'
numOrientations = 8
orientations = numpy.arange(0, 360,
                            360.0 / numOrientations)  # Remember, ranges in Python do NOT include the final value!
wedgeEdges = numpy.arange(0, 360+360.0/numWedges, 360.0/numWedges)
print wedgeEdges
startingPhase = 0  # initial phase for gratingStim

# create grating stims
gratingStim1 = visual.GratingStim(win=myWin, tex=textureType, units='deg', mask='circle', contrast=1,
                                  pos=centerPoint, size=stimSize, sf=spatialFreq, autoLog=False)
wedge1 = visual.RadialStim(myWin, units='deg', size=100, visibleWedge=[0,0], mask = [1])
wedge2 = visual.RadialStim(myWin, units='deg', size=100, visibleWedge=[0,0], mask = [1])


gratingStim1.setAutoDraw(True)
wedge1.setAutoDraw(True)
wedge2.setAutoDraw(True)
wedge1.setContrast(0)
wedge2.setContrast(0)

# Set up the trigger behavior
if triggerType == "NoTrigger":
    import noTrigger
    trigger = noTrigger.noTrigger(None)
elif triggerType == "SerialDaqOut" or triggerType == 'OutOnly':
    import serialTriggerDaqOut
    print 'Imported trigger serialTriggerDaqOut'
    trigger = serialTriggerDaqOut.serialTriggerDaqOut(serialPortName)
    # determine the Next experiment file name
    expName=trigger.getNextExpName([dataPath,date])
    print "Trial name: ",expName
    if triggerType == 'OutOnly':
        trigger.readSer=False
    if adjustDurationToMatch2P:
        print "Waiting for serial Triggers"
        stimDuration = trigger.extendStimDurationToFrameEnd(stimDuration)
    # store the stimulus data and prepare the directory
    trigger.preTrialLogging([dataPath,date,expName,stimCodeName,stiminfo,logFilePath])
elif triggerType=="DaqIntrinsicTrigger":
    import daqIntrinsicTrigger
    trigger = daqIntrinsicTrigger.daqIntrinsicTrigger(None)
else:
    print "Unknown trigger type", triggerType


# run
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
for trial in range(0, numTrials):
    stimWedgesOrder = range(0, numWedges+1)
    random.shuffle(stimWedgesOrder)
    # determine stim order
    print "Beginning Trial", trial + 1
    stimOrder = range(0, len(orientations))
    if isRandom:
        random.shuffle(stimOrder)
    for stimWedgesNumber in stimWedgesOrder:
        trigger.preStim(stimWedgesNumber + 1)
        if stimWedgesNumber == numWedges:
            gratingStim1.setContrast(0)
            print "\tStim", int(stimWedgesNumber + 1), 'blank'
        else:
            edge1 = wedgeEdges[stimWedgesNumber]
            print edge1
            edge2 = wedgeEdges[stimWedgesNumber+1]
            print edge2
            if edge1 < 179:
                if edge2 < 181:
                    wedge1.visibleWedge = [0, edge1+1]
                    wedge2.visibleWedge = [edge2, 360]
                else:
                    wedge1.visibleWedge = [0, 182]
                    wedge2.visibleWedge = [180, edge1]
            else:
                if edge2 > 359:
                    wedge1.visibleWedge = [0, 182]
                    wedge2.visibleWedge = [180, edge1]
                else:
                    wedge1.visibleWedge = [0, edge1]
                    wedge2.visibleWedge = [edge2, 360]
            print "\tStim", int(stimWedgesNumber + 1)
            gratingStim1.setContrast(contrast)
        for stimNumber in stimOrder:
            # convert orientations to standard lab notation
            gratingStim1.ori = orientations[stimNumber] - 90
            #print "\t", "\t", orientations[stimNumber], 'deg'
            clock.reset()
            while clock.getTime() < stimDuration:
                gratingStim1.setPhase(startingPhase + clock.getTime() * temporalFreq)
                trigger.preFlip(None)
                myWin.flip()
                trigger.postFlip(None)

        # now do ISI
        clock.reset()
        gratingStim1.setContrast(0)
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)
        while clock.getTime() < isi:
            # Keep flipping during the ISI. If you don't do this you can get weird flash artifacts when you resume flipping later.
            myWin.flip()
        trigger.postStim(None)

# trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'