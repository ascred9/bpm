# bpm
Beam profile monitor at g-2EDM
This macross build some hists and draw photos.
Short descriptions:
  read(filename, &pars) --- read the file and get parameters like a exposure time.
  read_and_fill(filename, factor) --- read and fill a 2d hist (photo), factor is used for recalculation pixels to mm. Return a TH2I hist pointer.
  build_1dimhist(hist2d) --- build a 1d hist with amplitudes of each pixel.
  calc_stat(array) --- calculate the mean and stddev of a dataset (array). It is used in other methods.
  fill_data(fname, &time, subtract_ped = true) --- read files with signal and pedestals, then subtract the pedestals. Return vector of amplitudes.
  fill_tree(dirname, ext, subtract_ped = true) --- read all files in the directory and fill tree with mean, stddev and exposure time.
  average(input, xsize, ysize) and filter(input, xsize, ysize) methods are used for filtering the input image to detect the ellipse edge of the target.
  fcn(npar, gin, f, par, iflag) --- function to fit the ellipse edge, is used by MINUIT. Calculate the distance between ellipse and pixela and summarize all of them.
  DrawEdges(filename) --- read a file, filter the input image to emphasize the ellipse edge, fit the pixels to get the ellipse parameters, and finnaly return the factor number to convert the pixels to mm's.
  fline(x, par) --- fit function of a straight line.
  DrawCanvas(filename, factor) --- draw the photo and also fit the profiles at the center of the picture. The factor number can be get from the above method (DrawEdges) or by the geometrical calculations from the distances. In this momemnt it doesn't subtract the pedestals.

The most usable case is to run the method DrawCanvas.
