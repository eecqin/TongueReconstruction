XRMBtools

(C) 2008 Chao Qin and Miguel A. Carreira-Perpinan
    Electrical Engineering and Computer Science
    University of California, Merced
    http://eecs.ucmerced.edu


The following are some Matlab functions useful to work with the Wisconsin
X-ray Microbeam database. For instructions of use see the individual files
and the summary below.

The Wisconsin X-ray Microbeam (XRMB) database is available from the
University of Wisconsin-Madison at http://www.medsch.wisc.edu/ubeam.


- XRMBdemo: example of use of the functions.

- XRMBread: reads X-ray microbeam data, acoustic wave and phonetic labels
  from an utterance file.

- XRMBscat: scatterplot of a set of X-ray frames in the midsaggital plane.

- XRMBtrace: displays traces of articulators for a given utterance. 

- XRMBframe: displays the 8 articulatory pellets and the interpolated
  tongue in the midsaggital plane for a single frame.

- XRMBdownsamp: performs downsampling on the original X-ray data.

- XRMBplot: displays selected acoustic features (e.g., acoustic wave,
  phonetic label, spectrogram, energy contour, pitch contour) and X-ray
  channels as a function of time, and optionally creates an animation.

- pcorrel: Pearson correlation between estimated and original (articulatory)
  trajectories.

- rmse: RMS error between estimated and original (articulatory) trajectories.


Some of these functions may require external functions, as follows:
- XRMBdownsamp needs 'resample' and 'filtfilt' from the Matlab Signal
  Processing Toolbox.
- XRMBplot needs 'hamming' from the Matlab Signal Processing Toolbox;
  Malcolm Slaney's 'MakeQTMovie' and 'ReadQTMovie'; and some functions from
  our speech_analysis toolbox, which can be obtained from the same web site
  as XRMBtools.

