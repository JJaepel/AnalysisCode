from psychopy import visual, logging, core, filters, event, monitors
import math, random, numpy, time, imp, sys

sys.path.append("Z:/Juliane/Code/Stimuli/triggers")  # path to trigger classes
from os import path

print "initialized"

# ---------- Stimulus Description ---------- #
'''A fullscreen drifting grating for 2pt orientation tuning'''
# ---------- Monitor Properties ----------#
mon = monitors.Monitor('testMonitor')  # gets the calibration
mon.setDistance(25)
overwriteGammaCalibration = False
newGamma = 0.479

# ---------- Stimulus Parameters ---------- #
# trials and duration
numTrials = 10  # Run all the stims this many times
doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.
nBlank = 1  # number of blanks to show per trial.
blankpercent = 0.1

numOrientations = 8
orientations = numpy.arange(0, 360,
                            360.0 / numOrientations)

textureType = 'sqr'  # 'sqr' = square wave, 'sin' = sinusoidal,'sqrDutyCycle'
if textureType == 'sqrDutyCycle':
    dutyCycle = 10  # can be 1, 2, 3, 4, 6, 8,
    textureType = 1 * numpy.ones((dutyCycle, 1));
    textureType[1, :] = -1;

stimDuration = .5
isi = 1

isRandom = 1
initialDelay = 1  # time in seconds to wait before first stimuli. Set to 0 to begin ASAP.

centerPoint = [0, 0]
stimSize = [5, 5]
stimPos = centerPoint
numStimElev = 4
numStimAzi = 8

startingPhase = 0.0
temporalFreq = 2.5
spatialFreq = 0.1

myMask = numpy.array([
        [ 1,1,1,1],
        [ 1,1,1,1],
        [1,1,1,1],
        [1,1,1,1],
        ])

barColor = 1 #1 for white, 0 for black, 0.5 for gray, 0 for black, etc.

if blankpercent > 0:
    nBlank = int(math.floor(numStimAzi * numStimElev * blankpercent))

# Triggering type
# Can be any of:#"NoTrigger" - no triggering; stim will run freely
# "SerialDaqOut" - Triggering by serial port. Stim codes are written to the MCC DAQ.
# "OutOnly" - no input trigger, but does all output (to CED) and logging
# "DaqIntrinsicTrigger" - waits for stimcodes on the MCC DAQ and displays the appropriate stim ID
# triggerType = 'DaqIntrinsicTrigger'
triggerType = 'NoTrigger'
serialPortName = 'COM2'  # ignored if triggerType is "None"
adjustDurationToMatch2P = True

# Experiment logging parameters
import time

dataPath = '//mpfiwdtfitz68/Spike2Data/'
date = (time.strftime("%Y-%m-%d"))
logFilePath = dataPath + date + '\\' + date + '.txt'  # including filepath

# ---------- Stimulus code begins here ---------- #
stimCodeName = path.dirname(path.realpath(__file__)) + '\\' + path.basename(__file__)

# make a window
myWin = visual.Window(size=[1920, 1080], monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=False)
if overwriteGammaCalibration:
    myWin.setGamma(newGamma)
    print "Overwriting Gamma Calibration. New Gamma value:", newGamma

print "made window, setting up triggers"
# Set up the trigger behavior
trigger = None
if triggerType == "NoTrigger":
    import noTrigger

    trigger = noTrigger.noTrigger(None)
elif triggerType == "SerialDaqOut" or triggerType == 'OutOnly':
    import serialTriggerDaqOut

    print 'Imported trigger serialTriggerDaqOut'
    trigger = serialTriggerDaqOut.serialTriggerDaqOut(serialPortName)
    # determine the Next experiment file name
    expName = trigger.getNextExpName([dataPath, date])
    print "Trial name: ", expName
    if triggerType == 'OutOnly':
        trigger.readSer = False
    # Record a bunch of serial triggers and fit the stim duration to an exact multiple of the trigger time
    if adjustDurationToMatch2P:
        print "Waiting for serial Triggers"
        stimDuration = trigger.extendStimDurationToFrameEnd(stimDuration)
    # store the stimulus data and prepare the directory
    trigger.preTrialLogging([dataPath, date, expName, stimCodeName])
elif triggerType == "DaqIntrinsicTrigger":
    import daqIntrinsicTrigger

    trigger = daqIntrinsicTrigger.daqIntrinsicTrigger(None)
else:
    print "Unknown trigger type", triggerType


gratingStim = visual.GratingStim(win=myWin, tex=textureType, units='deg',
                                 pos=centerPoint, size=stimSize, mask='circle', sf=spatialFreq, autoLog=False)
gratingStim.setAutoDraw(True)
barTexture = numpy.ones([256, 256, 3]);
minAzim = centerPoint[0] - ((numStimAzi * stimSize[0])/2) + stimSize[0]/2
minElev = centerPoint[1] + ((numStimElev * stimSize[1])/2) - stimSize[1]/2
print "made stims"
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
print "\n", str(numStimElev * numStimAzi + doBlank), "stims will be run for", str(numTrials), "trials."
if nBlank > 1:
    print "Will run blank " + str(nBlank) + " times"

duration = initialDelay + numTrials * numStimElev * numStimAzi * (isi + numOrientations*stimDuration) + numTrials * nBlank
print"estimated duration: " + str(math.ceil(duration/60)) + " min/" + str(duration) + " seconds"
if initialDelay > 0:
    print" waiting " + str(initialDelay) + " seconds before starting stim to acquire a baseline."
    gratingStim.setContrast(0)
    while clock.getTime() < initialDelay:
        myWin.flip()

for trial in range(0, numTrials):
    # determine stim order
    print "Beginning Trial", trial + 1
    stimOrder = range(0, numStimAzi * numStimElev + doBlank)
    if nBlank > 1:
        blankID = numStimAzi * numStimElev
        for ibl in range(1, nBlank):
            stimOrder.append(blankID)
    if isRandom:
        random.shuffle(stimOrder)
    #run

    for stimNumber in stimOrder:
        trigger.preStim(stimNumber + 1)
        if stimNumber == numStimAzi * numStimElev:
            print "\tStim", stimNumber + 1, " (blank)"
        else:
            xPos = minAzim + (stimNumber % numStimAzi) * stimSize[0]
            yPos = minElev - (math.ceil(stimNumber /numStimAzi)) * stimSize[1]
            print "\tStim", stimNumber + 1, ',', xPos, 'deg Azimuth,', yPos, 'deg Elevation'
        clock.reset()

        stimOris = range(0, numOrientations)
        random.shuffle(stimOris)
        for oris in stimOris:
            gratingStim.ori = orientations[oris] - 90
            while clock.getTime() < stimDuration:
                if stimNumber == numStimAzi * numStimElev:
                    gratingStim.setContrast(0)
                else:
                    gratingStim.setContrast(1)
                    gratingStim.pos = [xPos, yPos]
                    gratingStim.setPhase(startingPhase + clock.getTime() * temporalFreq)
                    trigger.preFlip(None)
                    myWin.flip()
                    trigger.postFlip(None)
            clock.reset()
        # now do ISI
        if isi != 0:
            clock.reset()
            print 'ISI'
            while clock.getTime() < isi:
                gratingStim.setContrast(0)
                trigger.preFlip(None)
                myWin.flip()
                trigger.postFlip(None)
        trigger.postStim(None)

#trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'