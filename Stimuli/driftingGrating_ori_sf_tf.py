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
numOrientations = 2  # typically 4, 8, or 16#
orientations = numpy.arange(0.0,180,180.0/numOrientations) #Remember, ranges in Python do NOT include the final value!

numTrials = 10  # Run all the stims this many times

doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.
nBlank = 1  # number of blanks to show per trial.
blankpercent = 0.1
changeDirectionAt = 1  # When do we change movement directions? If 1, there should be no reversal. Setting to 0 causes instantaneous reversal, effectively inverting all stimIDs and angles.
# DO NOT SET TO 0!.

stimDuration = 3
isi = 3

isRandom = 1
initialDelay = 5  # time in seconds to wait before first stimuli. Set to 0 to begin ASAP.

# Grating parameter
temporalFreq =(0.2, 0.8, 3.2, 6.4) #12.8 additional, 1.6 additional
spatialFreq = (0.04, 0.08, 0.16, 0.32) #0.02 and 0.06 additional
if blankpercent > 0:
    nBlank = int(math.floor(len(orientations) * len(spatialFreq) * len(temporalFreq) * blankpercent))
contrast = 1
textureType = 'sqr'  # 'sqr' = square wave, 'sin' = sinusoidal,'sqrDutyCycle'
if textureType == 'sqrDutyCycle':
    dutyCycle = 10  # can be 1, 2, 3, 4, 6, 8,
    textureType = 1 * numpy.ones((dutyCycle, 1))
    textureType[1, :] = -1
startingPhase = 0.0  # initial phase for gratingStim
# aperture and position parameters
centerPoint = [-65, 0]
stimSize = [500, 500]  # Size of grating in degrees

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
myWin = visual.Window(size=[3840, 1080], monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=False)
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
    import triggers.serialTriggerDaqOut

    print 'Imported trigger serialTriggerDaqOut'
    trigger = triggers.serialTriggerDaqOut.serialTriggerDaqOut(serialPortName)
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
    trigger.preTrialLogging([dataPath, date, expName, stimCodeName, orientations, logFilePath])
elif triggerType == "DaqIntrinsicTrigger":
    import triggers.daqIntrinsicTrigger

    trigger = triggers.daqIntrinsicTrigger.daqIntrinsicTrigger(None)
else:
    print "Unknown trigger type", triggerType

print stimDuration
changeDirectionTimeAt = stimDuration * changeDirectionAt

# create grating stim
gratingStim = visual.GratingStim(win=myWin, tex=textureType, units='deg',
                                 pos=centerPoint, size=stimSize, mask='none', sf=spatialFreq[1], autoLog=False)
gratingStim.setAutoDraw(True)
barTexture = numpy.ones([256, 256, 3])
# flipStim = visual.PatchStim(win=myWin,tex=barTexture,mask='none',units='pix',pos=[-920,500],size=(100,100))
# flipStim.setAutoDraw(True)#up left, this is pos in y, neg in x
clrctr = 1
print "made grating"
# run
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
print "\n", str(len(orientations) * len(spatialFreq) * len(temporalFreq) + doBlank), "stims will be run for", str(numTrials), "trials."
if nBlank > 1:
    print "Will run blank " + str(nBlank) + " times"
# force a wait period of at least 5 seconds before first stim
duration = initialDelay + numTrials * numOrientations * len(spatialFreq) * len(temporalFreq) * (
            isi + stimDuration) + numTrials * nBlank
print"estimated duration: " + str(math.ceil(duration / 60)) + " min/" + str(duration) + " seconds"
if initialDelay > 0:
    print" waiting " + str(initialDelay) + " seconds before starting stim to acquire a baseline."
    gratingStim.setContrast(0)
    # flipStim.setContrast(0)
    while clock.getTime() < initialDelay:
        myWin.flip()
for trial in range(0, numTrials):
    # determine stim order
    print "Beginning Trial", trial + 1
    stimOrder = range(0, len(orientations) * len(temporalFreq) * len(spatialFreq) + doBlank)
    if nBlank > 1:
        blankID = len(orientations)*len(spatialFreq)*len(temporalFreq)
        for ibl in range(1, nBlank):
            stimOrder.append(blankID)
    if isRandom:
        random.shuffle(stimOrder)
    for stimNumber in stimOrder:
        # display each stim
        trigger.preStim(stimNumber + 1)
        # display stim
        # flipStim.setContrast(1)
        # flipStim.setAutoDraw(True)
        # convert orientations to standard lab notation
        if stimNumber == len(orientations)*len(spatialFreq)*len(temporalFreq):
            gratingStim.setContrast(0)
            print "\tStim", stimNumber + 1, " (blank)"
        else:
            gratingStim.setContrast(contrast)
            sf = math.ceil(stimNumber / (len(temporalFreq) * len(orientations)))
            newNumber = stimNumber - (int(sf)) * len(temporalFreq) * len(orientations)
            ori = newNumber % len(orientations)
            tf = math.ceil(newNumber / len(orientations))
            gratingStim.ori = orientations[int(ori)] - 90
            gratingStim.setSF(spatialFreq[int(sf)])
            print "\tStim", stimNumber + 1, orientations[int(ori)], 'deg', temporalFreq[int(tf)], 'Hz', spatialFreq[int(sf)], 'cpd'
        clock.reset()
        while clock.getTime() < stimDuration:
            clrctr = clrctr + 1
            # if clrctr%2==1:
            # flipStim.setColor((0,0,0),colorSpace='rgb')
            # flipStim.setContrast(-1)
            # else:
            # flipStim.setColor((1,1,1),colorSpace='rgb')
            # flipStim.setContrast(1)
            if clock.getTime() > changeDirectionTimeAt:
                gratingStim.setPhase(startingPhase + changeDirectionTimeAt * temporalFreq[int(tf)] - (
                            clock.getTime() - changeDirectionTimeAt) * temporalFreq[int(tf)])
                # print startingPhase+changeDirectionTimeAt*temporalFreq - (clock.getTime()-changeDirectionTimeAt)*temporalFreq
            else:
                if stimNumber == len(orientations) * len(temporalFreq):
                    gratingStim.setPhase(startingPhase + clock.getTime() * temporalFreq[1])
                else:
                    gratingStim.setPhase(startingPhase + clock.getTime() * temporalFreq[int(tf)])
            trigger.preFlip(None)
            myWin.flip()
            trigger.postFlip(None)

        # now do ISI
        if isi != 0:
            clock.reset()
            print 'ISI'
            gratingStim.setContrast(0)
            # flipStim.setContrast(0)
            trigger.preFlip(None)
            myWin.flip()
            trigger.postFlip(None)
            # flipStim.setAutoDraw(False)
            while clock.getTime() < isi:
                # Keep flipping during the ISI. If you don't do this you can get weird flash artifacts when you resume flipping later.
                myWin.flip()
        trigger.postStim(None)

#trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'