MOCHAtools

(C) 2008 Chao Qin and Miguel A. Carreira-Perpinan
    Electrical Engineering and Computer Science
    University of California, Merced
    http://eecs.ucmerced.edu


The following are some Matlab functions useful to work with the
MOCHA-TIMIT database. For instructions of use see the individual files and
the summary below.

The MOCHA-TIMIT database is available from the Centre for Speech
Technology Research (CSTR) of the University of Edinburgh at
http://www.cstr.ed.ac.uk/research/projects/artic/mocha.html.


- MOCHAdemo: example of use of the functions.

- MOCHAread: reads EMA data, acoustic wave and phonetic labels from an
  utterance file.

- MOCHAscat: scatterplot of a set of EMA frames in the midsaggital plane.

- MOCHAtrace: displays traces of articulators for a given utterance. 

- MOCHAframe: displays the articulatory pellets and the interpolated
  tongue in the midsaggital plane for a single frame.

- MOCHAdownsamp: performs downsampling on the original EMA data.

- MOCHAplot: displays selected acoustic features (e.g., acoustic wave,
  phonetic label, spectrogram, energy contour, pitch contour) and EMA
  channels as a function of time, and optionally creates an animation.

- pcorrel: Pearson correlation between estimated and original (articulatory)
  trajectories.

- rmse: RMS error between estimated and original (articulatory) trajectories.


Some of these functions may require external functions, as follows:
- MOCHAdownsamp needs 'resample' and 'filtfilt' from the Matlab Signal
  Processing Toolbox.
- MOCHAplot needs 'hamming' from the Matlab Signal Processing Toolbox;
  Malcolm Slaney's 'MakeQTMovie' and 'ReadQTMovie'; and some functions from
  our speech_analysis toolbox, which can be obtained from the same web site
  as MOCHAtools.

