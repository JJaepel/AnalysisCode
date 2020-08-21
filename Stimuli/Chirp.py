from psychopy import visual, logging, core, filters, event, monitors
import math, random, numpy, time, imp, sys
import matplotlib.pyplot as plt

sys.path.append("../triggers") #path to trigger classes
from os import path

print "initialized"

#Experiment logging parameters
dataPath='X:/'
animalName='F871_2018-06-21'; 

# ---------- Stimulus Description ---------- #
'''A fullscreen chirp stimulus for dLGN'''
# ---------- Monitor Properties ----------#
mon = monitors.Monitor('testMonitor')  # gets the calibration
mon.setDistance(25)
overwriteGammaCalibration = False
newGamma = 0.479

numTrials = 10 # ca. 30 s per Trial
initialDelay = 5
stimDuration = 3
stimDuration_short = 2
isi = 2

chirpDur = 6
chirpMaxFreq = 18
ContrastFreq = 5.4
durFr = 1/80.0

chirps = int(chirpDur / durFr)
# print"Chirps: " + str(chirps)
K_HzPerSec = chirpMaxFreq / chirpDur

centerPoint = [0, 0]
stimSize = [500, 500]

# Triggering type
# Can be any of:#"NoTrigger" - no triggering; stim will run freely
# "SerialDaqOut" - Triggering by serial port. Stim codes are written to the MCC DAQ.
# "OutOnly" - no input trigger, but does all output (to CED) and logging
# "DaqIntrinsicTrigger" - waits for stimcodes on the MCC DAQ and displays the appropriate stim ID
# triggerType = 'DaqIntrinsicTrigger'
triggerType = 'OutOnly'
serialPortName = 'COM2'  # ignored if triggerType is "None"
adjustDurationToMatch2P = True

date = (time.strftime("%Y-%m-%d"))
logFilePath = dataPath + date + '\\' + date + '.txt'  # including filepath

# ---------- Stimulus code begins here ---------- #
stimCodeName = path.dirname(path.realpath(__file__)) + '\\' + path.basename(__file__)

myWin = visual.Window(size=[1920, 1080], monitor=mon, fullscr=True, screen=1, allowGUI=False, waitBlanking=False)
print "made window, setting up triggers"
#Set up the trigger behavior
trigger = None
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
        flashInterval = trigger.extendStimDurationToFrameEnd(flashInterval)
    # store the stimulus data and prepare the directory
    trigger.preTrialLogging([dataPath,animalName,expName,stimCodeName,orientations,logFilePath])
elif triggerType=="DaqIntrinsicTrigger":
    import daqIntrinsicTrigger
    trigger = daqIntrinsicTrigger.daqIntrinsicTrigger(None) 
else:
    print "Unknown trigger type", triggerType
        

# create fullcontrast and halfcontrast
Halfval = 127
Dark = visual.PatchStim(myWin, tex="sqrXsqr",texRes=64, color= 0, colorSpace= 'rgb255',
    size=stimSize, sf=.0008, mask = 'none', pos=centerPoint)
Full = visual.PatchStim(myWin, tex="sqrXsqr",texRes=64, color= 254, colorSpace= 'rgb255',
    size=stimSize, sf=.0008, mask = 'none', pos=centerPoint)
Half = visual.PatchStim(myWin, tex="sqrXsqr",texRes=64, color= Halfval , colorSpace= 'rgb255',
    size=stimSize, sf=.0008, mask = 'none', pos=centerPoint)

print "made stims"
# run
clock = core.Clock()  # make one clock, instead of a new instance every time. Use
duration = initialDelay + numTrials * (3*stimDuration_short+ 2*stimDuration+ 2*(chirpDur+2)+ isi)
print"estimated duration: " + str(math.ceil(duration / 60)) + " min/" + str(duration) + " seconds"
if initialDelay > 0:
    print" waiting " + str(initialDelay) + " seconds before starting stim to acquire a baseline."
    Half.draw()
    while clock.getTime() < initialDelay:
        myWin.flip()

for trial in range(0,numTrials):
    print "Beginning Trial", trial + 1
    # ON and OFF stimuli
    clock.reset()
    while clock.getTime() < stimDuration_short:
        Dark.draw()
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)
    clock.reset()

    print 'ON'
    trigger.preStim(1)
    while clock.getTime() < stimDuration:
        Full.draw()
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)
    clock.reset()

    print 'OFF'
    trigger.preStim(2)
    while clock.getTime() < stimDuration:
        Dark.draw()
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)
    clock.reset()

    trigger.preStim(0)
    while clock.getTime() < stimDuration_short:
        Half.draw()
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)

    # Frequency chirp
    print 'frequency chirp'
    trigger.preStim(3)
    for freq in range(0, chirps):
        t_s = freq * durFr
        Intensity = math.sin(math.pi * K_HzPerSec * t_s**2) * Halfval + Halfval
        clock.reset()
        while clock.getTime() < durFr:
            Int = visual.PatchStim(myWin, tex="sqrXsqr", texRes=64, color=Intensity, colorSpace='rgb255',
                                    size=stimSize, sf=.0008, mask='none', pos=centerPoint)
            Int.draw()
            trigger.preFlip(None)
            myWin.flip()
            trigger.postFlip(None)
    clock.reset()
    # Gap between frequency chirp and contrast chirp

    trigger.preStim(0)
    while clock.getTime() < stimDuration_short:
        Half.draw()
        trigger.preFlip(None)
        myWin.flip()
        trigger.postFlip(None)

    # Contrast chirp
    print 'contrast chirp'
    trigger.preStim(4)
    for cont in range(0, chirps):
        t_s = cont * durFr
        IRamp = int(Halfval * t_s / chirpDur)
        IntensityCon = math.sin(math.pi * ContrastFreq * t_s ) * IRamp + Halfval
        clock.reset()
        while clock.getTime() < durFr:
            Int = visual.PatchStim(myWin, tex="sqrXsqr", texRes=64, color=IntensityCon, colorSpace='rgb255',
                                   size=stimSize, sf=.0008, mask='none', pos=centerPoint)
            Int.draw()
            trigger.preFlip(None)
            myWin.flip()
            trigger.postFlip(None)

    # Gap after contrast chirp
    trigger.preStim(0)
    if isi != 0:
        clock.reset()
        print 'ISI'
        while clock.getTime() < isi:
            Half.draw()
            trigger.preFlip(None)
            myWin.flip()
            trigger.postFlip(None)

    trigger.postStim(None)

#trigger.wrapUp([logFilePath, expName])
print 'Finished all stimuli.'