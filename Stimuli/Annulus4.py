from psychopy import visual, logging, core, filters, event, monitors, tools
import pylab, math, random, numpy, time, imp, sys, itertools

sys.path.append("Z:/Juliane/Code/Stimuli/triggers")  # path to trigger classes
from os import path

print "initialized"

# ---------- Stimulus Description ---------- #
'''A fullscreen drifting grating for 2pt orientation tuning'''
# ---------- Monitor Properties ----------#
mon = monitors.Monitor('testMonitor')  # gets the calibration for stimMonitor
mon.setDistance(25)
myWin = visual.Window(size=mon.getSizePix(), monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=True,
                      units='deg', blendMode='add')

# aperture and position parameters
# centerpos = tools.monitorunittools.pix2deg( myWin.size[0]/4,mon)#20#68 #visual degrees off by 2 (33.15 x 2)\
# centerPoint = [-centerpos,0] #[0,0]
centerPoint = [0, 0]
stimSize1 = [66, 66]  # [300,300]Size of grating in degrees
stimSize2 = [22, 22]

# ---------- Stimulus Parameters ---------- #
# trials and duration

doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.

numOrientations = 8  # typically 4, 8, or 16
orientations = numpy.arange(0.0, 360,
                            float(360 / numOrientations))  # Remember, ranges in Python do NOT include the final value!
orientations2 = numpy.arange(0.0, 360,
                             float(360 / numOrientations))  # Remember, ranges in Python do NOT include the final value!

allStimPairs = list(itertools.product(orientations, orientations2))  # combine stim conditions

numTrials = 6  # Run all the stims this many times
pauseBetweenTrials = 0  # Wait for user to press space bar in between trials
doBlank = 1  # 0 for no blank stim, 1 to have a blank stim. The blank will have the highest stimcode.
stimDuration = 2  # stimulus duration in seconds; will be adjusted if adjustDurationToMatch2P=1
isMoving = 1  # 1 for moving stimulus, 0 for stationary
changeDirectionAt = 1  # When do we change movement directions? If 1, there should be no reversal.  If 0.5, then movement reverses directions at stimDuration/2
isi = 4  # Kuo 4 to 2
isRandom = 1
# Grating parameter
temporalFreq = 4
spatialFreq = 0.2
contrast = 1  # 0.125
textureType = 'sqr'  # 'sqr' = square wave, 'sin' = sinusoidal
startingPhase = 0  # initial phase for gratingStim
changeDirectionAt = 1;  # 0.5#When do we change movement directions? If 1, there should be no reversal. Setting to 0 causes instantaneous reversal, effectively inverting all stimIDs and angles.
changeDirectionTimeAt = stimDuration * changeDirectionAt

# create grating stim
gratingStim1 = visual.GratingStim(win=myWin, tex=textureType, units='deg', mask='circle', contrast=1,
                                  pos=centerPoint, size=stimSize1, sf=spatialFreq, autoLog=False)
gratingStim2 = visual.GratingStim(win=myWin, tex=textureType, units='deg', mask='circle', opacity=1, contrast=1,
                                  pos=centerPoint, size=stimSize2, sf=spatialFreq, autoLog=False)

gratingStim1.setAutoDraw(True)
gratingStim2.setAutoDraw(True)

barTexture = numpy.ones([256, 256, 3]);
flipStim = visual.PatchStim(win=myWin, tex=barTexture, mask='none', units='pix', pos=[-920, 500], size=(100, 100))
flipStim.setAutoDraw(True)  # up left, this is pos in y, neg in x
clrctr = 1;

triggerType = 'NoTrigger'
serialPortName = 'COM2'  # ignored if triggerType is "None"
adjustDurationToMatch2P = True


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
    trigger.preTrialLogging([dataPath, date, expName, stimCodeName, orientations, logFilePath])
elif triggerType == "DaqIntrinsicTrigger":
    import daqIntrinsicTrigger

    trigger = daqIntrinsicTrigger.daqIntrinsicTrigger(None)
else:
    print "Unknown trigger type", triggerType
print stimDuration

# run the stimuli in random order
clrctr = 1;
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
print "\n", str(len(allStimPairs) + doBlank), "stims will be run for", str(numTrials), "trials."
print "\nTotal duration will be", (1 / 60.0) * (len(allStimPairs) + doBlank) * numTrials * (
            stimDuration + isi), "minutes."

# run
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
for trial in range(0, numTrials):
    # determine stim order
    if triggerType != "NoTrigger":
        if pauseBetweenTrials:
            print "Waiting for User Input"
            event.waitKeys(keyList=['space'])
    print "Beginning Trial", trial + 1
    stimOrder = range(0, len(allStimPairs) + doBlank)
    if isRandom:
        random.shuffle(stimOrder)
    for stimNumber in stimOrder:

        # display each stim
        trigger.preStim(stimNumber + 1)
        # display stim
        flipStim.setContrast(1)
        flipStim.setAutoDraw(True)
        # convert orientations to standard lab notation
        if stimNumber == len(allStimPairs):
            gratingStim1.setContrast(0)
            gratingStim2.setContrast(0)
            print "\tStim", stimNumber + 1, " (blank)"
        else:
            gratingStim1.setContrast(contrast)
            gratingStim2.setContrast(0)
            gratingStim1.ori = allStimPairs[stimNumber][0] - 90.00
            gratingStim2.ori = allStimPairs[stimNumber][1] - 90.00
            print "\tStim", stimNumber + 1, allStimPairs[stimNumber], 'deg'
        clock.reset()
        while clock.getTime() < stimDuration:
            clrctr = clrctr + 1;
            if clrctr % 2 == 1:
                # flipStim.setColor((0,0,0),colorSpace='rgb')
                flipStim.setContrast(-1)
            else:
                # flipStim.setColor((1,1,1),colorSpace='rgb')
                flipStim.setContrast(1)
            if isMoving:
                if clock.getTime() > changeDirectionTimeAt:
                    gratingStim1.setPhase(startingPhase + changeDirectionTimeAt * temporalFreq - (
                                clock.getTime() - changeDirectionTimeAt) * temporalFreq)
                    gratingStim2.setPhase(startingPhase + changeDirectionTimeAt * temporalFreq - (
                                clock.getTime() - changeDirectionTimeAt) * temporalFreq)
                else:
                    gratingStim1.setPhase(startingPhase + clock.getTime() * temporalFreq)
                    gratingStim2.setPhase(startingPhase + clock.getTime() * temporalFreq)
            trigger.preFlip(None)
            myWin.flip()
            trigger.postFlip(None)

        # now do ISI
        clock.reset()
        gratingStim1.setContrast(0)
        gratingStim2.setContrast(0)
        flipStim.setContrast(0)
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)
        flipStim.setAutoDraw(False)
        while clock.getTime() < isi:
            # Keep flipping during the ISI. If you don't do this you can get weird flash artifacts when you resume flipping later.
            myWin.flip()
        trigger.postStim(None)

# trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'