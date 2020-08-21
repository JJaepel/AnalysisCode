import serial, csv, time, math, datetime
from psychopy import visual, core, filters, monitors
import pylab, math, random, numpy, sys

sys.path.append("Z:/Juliane/Code/Stimuli/triggers")
#import UniversalLibrary as UL
from abstractTrigger import trigger
from os import path,makedirs
import shutil
import glob



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

mon = monitors.Monitor('testMonitor')  # gets the calibration

# ---------- Stimulus Parameters ---------- #
numOrientations = 16
orientations = numpy.arange(0, 360,
                            360.0 / numOrientations)  # Remember, ranges in Python do NOT include the final value!
stiminfo = orientations

numTrials = 10

doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.
changeDirectionAt = 0  # When do we change movement directions? If 1, there should be no reversal. Setting to 0 causes instantaneous reversal, effectively inverting all stimIDs and angles.
stimDuration = 0.15
isi = 0.5
isRandom = 1

# 24 deg increments, first stimulus no annulus, rest use a 50% annulus

stimSizes = [5, 10, 20, 40, 60, 80]

stimPos = [0, 0]

temporalFreq = 4
spatialFreq = 0.1
contrast = 1

textureType = 'sqr'  # 'sqr' = square wave, 'sin' = sinusoidal, 'sqrDutyCycle'

# ---------- Stimulus code begins here ---------- #
stimCodeName = path.dirname(path.realpath(__file__)) + '\\' + path.basename(__file__)

# make a window
myWin = visual.Window(size=[1920, 1080], monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=False)
print "made window, setting up triggers"

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

print stimDuration

radius = filters.makeRadialMatrix(256)
# note - creating some space in between annuli (its not 0.5)
annulus = numpy.where(radius < 1, 1, 0) * numpy.where(radius > 0.5, 1, 0) * 2 - 1

# create grating stim
gratingStim = visual.GratingStim(win=myWin, mask=annulus, tex=textureType, units='deg',
                                 pos=stimPos, sf=spatialFreq, autoLog=False)
gratingStim.setAutoDraw(True)

barTexture = numpy.ones([256, 256, 3]);
flipStim = visual.PatchStim(win=myWin, tex=barTexture, mask='None', units='pix', pos=[-1920, 540], size=(100, 100))
flipStim.setAutoDraw(True)  # up left, this is pos in y, neg in x
clrctr = 1;
print "made grating"

stimOrder = range(0, len(orientations))

stimAnnulusOrder = range(0, len(stimSizes))

# run
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
print "\n", len(stimOrder), "stims will be run for", str(numTrials), "trials."

for trial in range(0, numTrials):

    # determine stim order
    print "Beginning Trial", trial + 1
    if isRandom:
        random.shuffle(stimOrder)

    for block in range(0, 3):

        if block == 0:
            print "center stim"
            gratingStim.mask = 'circle';
            stimbuffer = 0;
        elif block == 1:
            print "annulus (50%)"
            stimbuffer = 5;  # number of sizes used
            gratingStim.mask = annulus;
        elif block == 2:
            # blank
            stimbuffer = 5 * 2;

        while clock.getTime() < 2:
            myWin.flip()

        if block == 2:
            trigger.preStim(1 + stimbuffer)
            print "\tStim", 1 + stimbuffer, " (blank)"
            gratingStim.setContrast(0)
            clock.reset()
            while clock.getTime() < stimDuration:
                clrctr = clrctr + 1;
                if clrctr % 2 == 1:
                    flipStim.setContrast(-1)
                else:
                    flipStim.setContrast(1)
                gratingStim.setPhase(clock.getTime() * temporalFreq)
                trigger.preFlip(None)
                myWin.flip()
                trigger.postFlip(None)
        else:

            for stimAnnulusNumber in stimAnnulusOrder:

                trigger.preStim(stimAnnulusNumber + 1 + stimbuffer)

                gratingStim.size = [stimSizes[stimAnnulusNumber], stimSizes[stimAnnulusNumber]]

                print "\tStim", stimAnnulusNumber + 1 + stimbuffer, stimSizes[stimAnnulusNumber], 'deg'

                for stimNumber in stimOrder:

                    # display stim
                    flipStim.setContrast(1)
                    flipStim.setAutoDraw(True)
                    gratingStim.setContrast(contrast)
                    gratingStim.ori = orientations[stimNumber] - 90
                    clock.reset()
                    while clock.getTime() < stimDuration:
                        clrctr = clrctr + 1;
                        if clrctr % 2 == 1:
                            flipStim.setContrast(-1)
                        else:
                            flipStim.setContrast(1)
                        gratingStim.setPhase(clock.getTime() * temporalFreq)
                        trigger.preFlip(None)
                        myWin.flip()
                        trigger.postFlip(None)

                # now do ISI
                clock.reset()
                gratingStim.setContrast(0)
                flipStim.setContrast(0)
                trigger.preFlip(None)
                myWin.flip()
                trigger.postFlip(None)
                flipStim.setAutoDraw(False)
                while clock.getTime() < isi:
                    # Keep flipping during the ISI. If you don't do this you can get weird flash artifacts when you resume flipping later.
                    myWin.flip()
                trigger.postStim(None)

#trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'