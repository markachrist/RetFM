import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Formatter;
import java.util.Hashtable;
import java.util.Scanner;

import javax.swing.JFileChooser;

import ij.IJ;
import ij.ImagePlus;
import ij.Menus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.io.FileSaver;
import ij.io.OpenDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.EDM;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.frame.RoiManager;
import ij.process.BinaryProcessor;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

/**
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * Computes both shape and intensity statistics of retinal cells. The cells are first extracted using a thresholding and 
 * watershed segmentation approach. This approach also includes filters that exclude cells outside of a specified size range
 * and those outside circularity constraints. After segmentation, ROIs corresponding to the cells are added to the input image.
 * A set of shape and color/intensity features are computed for each cell and displayed in a results table. The set of 
 * measurements includes shape, intensity, and color features.
 * 
 * Please see the following publication for more details and cite the publication if using this plugin/code in any research.
 * 
 * Hedberg-Buenz et al. 2015. RetFM-J, an ImageJ-based module for automated counting and quantifying features of nuclei in 
 * retinal whole-mounts. Exp Eye Res.
 * 
 * Copyright Mark Christopher & The University of Iowa
 *
 */
public class RetFM_J implements PlugIn {
	
	/** Scale factor mapping pixels to length units */
	double scale = 3.15;
	
	/** Minimum size (area) of cells to count */
	double sizeMin = 10;
	
	/** Maximum size (area) of cells to count */
	double sizeMax = 150;
	
	/** Minimum circularity of cells to count */
	double circMin = 0.0;
	
	/** Maximum circularity of cells to count */
	double circMax = 1.0;
	
	/** Block size used for CLAHE local contrast adjustment */
	int contrastBlockSize = 64;
	
	/** Boolean indicating whether cells should be counted inside/outside the current ROI */
	boolean invertRoi = true;
	
	/** Boolean indicating whether CLAHE local contrast adjustment should be performed */
	boolean runContrast = true;
	
	/** Boolean indicating whether the plugin should display results for each image and wait for user input */
	boolean pauseBetweenImages = false;
	
	/** Boolean indicating whether the plugin should output to a single file or one file per image */
	boolean combineOutput = true;
	
	/** Boolean indicating whether the plugin should save images overlaying the detected cells */
	boolean saveOverlays = true;
	
	/** Boolean indicating whether the plugin should display total & examined image areas on completion */
	boolean displayAreas = true;
	
	/** Booleans indicating whether each channel should be considered */
	boolean useChannels[];
	
	/** The commands that can be called to run CLAHE contrast adjustment */
	String fijiContrast = "Enhance Local Contrast (CLAHE)";
	String ijContrast = "CLAHE ";
	
	/** The args that can be called to run CLAHE contrast adjustment */
	String fijiArgs = "blocksize=BLOCKSIZE histogram=256 maximum=3 mask=*None*";
	String ijArgs = "blocksize=BLOCKSIZE histogram=256 maximum=3";
	String contrastArgs = ijArgs;
	
	/** Set of feature image commands to compute for the images (requires FeatureJ plugin) */
	String featureCmds[] = {
			"FeatureJ Derivatives", "FeatureJ Derivatives", "FeatureJ Derivatives",
			"FeatureJ Edges", "FeatureJ Edges", "FeatureJ Edges",
			"FeatureJ Hessian", "FeatureJ Hessian", "FeatureJ Hessian",
			"FeatureJ Laplacian", "FeatureJ Laplacian", "FeatureJ Laplacian"
			};
	
	/** Arguments passed to compute the corresponding features indicated by featureCmds */
	String featureArgs[] = {
			"x-order=0 y-order=0 z-order=0 smoothing=1.0", 
			"x-order=0 y-order=0 z-order=0 smoothing=2.0", 
			"x-order=0 y-order=0 z-order=0 smoothing=4.0",
			"compute smoothing=1.0 lower=[] higher=[]", 
			"compute smoothing=2.0 lower=[] higher=[]", 
			"compute smoothing=4.0 lower=[] higher=[]",
			"largest absolute smoothing=4.0", 
			"largest absolute smoothing=4.0",
			"largest absolute smoothing=4.0", 
			"compute smoothing=1.0",
			"compute smoothing=2.0",
			"compute smoothing=4.0"
			};
	
	/** String used to delimit columns in ResultsTable row strings */
	String tableDelim = "\t";
	
	/** String used to delimit columns in output files */
	String outDelim = ",";
	
	/** Shape/location measurements to collect for each segmented cell */
	String measurements = "area centroid center perimeter bounding fit shape feret's area_fraction stack";
	
	/** Color measurements to collect for each segmented cell */
	String colorMeasurements = "mean standard modal min integrated median skewness kurtosis";
	
	/** Directory to store output */
	String outDir = null;
	
	/** Output buffer to store features computed from images */
	String outputCellFeatures = "";
	
	/**
	 * Performs the cell counting analysis on all currently open image windows.
	 * 
	 * See {publication} for the details of this analysis. 
	 *  
	 *  @param arg
	 *    Ignored argument string.
	 */
	public void run(String arg) {
		int ids[] = WindowManager.getIDList();
		
		if(this.createDialog()){
			
			boolean contrast = this.runContrast;
			double totalArea, examArea;
			String displayNames[] = new String[ids.length];
			String names[] = new String[ids.length];
			double totalAreas[] = new double[ids.length];
			double examAreas[] = new double[ids.length];
			
			//Constructor RoiManager(boolean hideWindow) never displays the window...
			RoiManager rois = new RoiManager(true);
			
			IJ.showProgress(0.0);

			//For every currently open image
			for(int i = 0; i < ids.length; ++i){
				
				int numCells = 0;
				String data[];
				
				//Get the current image/ROI
				ImagePlus p = WindowManager.getImage(ids[i]);
				ImagePlus color = p;
				
				//Create filename for text output
				displayNames[i] = color.getTitle();
				names[i] = color.getTitle() + ".txt";
				
				//Get current image and bring to front
				WindowManager.setCurrentWindow(p.getWindow());
				WindowManager.toFront(p.getWindow());
				
				String curTitle = p.getTitle();
				Roi r = p.getRoi();

				//Invert ROI, if needed
				if(invertRoi && r != null){
					ShapeRoi selection = new ShapeRoi(r);
					ShapeRoi all = new ShapeRoi(new Roi(0,0, p.getWidth(), p.getHeight()));
					r = all.xor(selection);
					WindowManager.getCurrentImage().setRoi(r);
				}
				
				//Run contrast adjustment
				if(contrast){
					contrast = this.adjustContrast(contrastBlockSize);
					
					if(!contrast){
						String m = "No contrast adjustment available.\n";
						m += "Install CLAHE (http://rsbweb.nih.gov/ij/plugins/clahe/index.html) ";
						m += "to include this functionality";
						IJ.showMessage("ImageJ", m);
					}
				}
				
				//Convert to gray
				FloatProcessor grayp = this.convertToGrayscale(p.getProcessor().convertToColorProcessor());
				ImageProcessor ip = grayp.convertToByte(false);
				ip.setRoi(r);
				
				//Find grayscale threshold via Otsu method
				int t = this.getThreshold(ip);
				ip.threshold(t);
				ip.invert();
				
				//Perform watershed segmentation
				EDM dm = new EDM();
				dm.toWatershed(ip);
				
				//Convert to binary
				BinaryProcessor bp = new BinaryProcessor((ByteProcessor)ip);
				bp.invert();
				p = new ImagePlus(curTitle, bp);
				p.setRoi(r);

				//Set pixel scaling based on user input
				Calibration calib = new Calibration();
				calib.pixelWidth = 1.0 / scale;
				calib.pixelHeight = 1.0 / scale;
				p.setCalibration(calib);
				
				//Get area info and write it out, if applicable
				totalArea = p.getWidth() * p.getHeight() / (scale * scale);
				examArea = p.getStatistics().area;
				
				totalAreas[i] = totalArea;
				examAreas[i] = examArea;
				
				//Analyze/count particles
				IJ.run("Set Measurements...", measurements + "redirect=None decimal=3");
				ResultsTable results = new ResultsTable();
				ParticleAnalyzer.setResultsTable(results);

				ParticleAnalyzer.setRoiManager(rois);
				
				//Set up particle analyzer and count cells
				int partOpts = ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES | ParticleAnalyzer.SHOW_NONE;
				ParticleAnalyzer analyze = new ParticleAnalyzer(partOpts, 0, null, 
						sizeMin*scale*scale, sizeMax*scale*scale,
						circMin, circMax);
				analyze.analyze(p);
				
				//Get shape/location measurements
				rois.runCommand("Show All with labels");
				rois.runCommand("Measure");

				//Copy data from results table
				numCells = results.getCounter();

				data = new String[numCells + 1];
				data[0] = results.getColumnHeadings();
				
				//Add Cell_ID column name
				data[0] = data[0].substring(1);
				data[0] += tableDelim;
				data[0] = "Cell_ID" + data[0];
				for(int row = 0; row < numCells; ++row){
					String rowStr = results.getRowAsString(row);
					data[row + 1] = rowStr;
				}

				ResultsTable.getResultsWindow().close(false);
				
				//Get measurements for each color channel
				IJ.run("Set Measurements...", colorMeasurements + "redirect=None decimal=3");
				
				for(int c = 0; c < this.useChannels.length; ++c){
					
					//Isolate current channel
					ByteProcessor channel = new ByteProcessor(color.getWidth(), color.getHeight(), 
							color.getProcessor().convertToColorProcessor().getChannel(c + 1));
					ImagePlus chanP = new ImagePlus("c_" + c, channel);
					
					WindowManager.setTempCurrentImage(chanP);
					
					//Take measurements
					rois.runCommand("Measure");
					results = ResultsTable.getResultsTable();
					
					//Get column names names
					String heads[] = results.getColumnHeadings().split(tableDelim);
					for(int h = 1; h < heads.length; ++h){
						data[0] += "c" + c + "_" + heads[h] + tableDelim;
					}
					
					//Get measurements from table
					for(int row = 0; row < numCells; ++row){
						String rowStr = results.getRowAsString(row);
						rowStr = rowStr.substring(rowStr.indexOf(tableDelim) + 1);
						data[row + 1] += tableDelim + rowStr;
					}
					
					ResultsTable.getResultsWindow().close(false);
				}
				
				/*******************************************************************************/
				//Get measurements for each feature image
				ImagePlus gray = new ImagePlus("gray", grayp);
				gray.show();
				WindowManager.setCurrentWindow(gray.getWindow());
				
				for(int f = 0; f < featureCmds.length; ++f){
					
					//Compute feature image
					IJ.run(featureCmds[f], featureArgs[f]);
					
					//Take measurements from features
					rois.runCommand("Measure");
					results = ResultsTable.getResultsTable();
					
					//Get column names names
					String heads[] = results.getColumnHeadings().split(tableDelim);
					for(int h = 1; h < heads.length; ++h){
						data[0] += "f" + f + "_" + heads[h] + tableDelim;
					}
					
					//Get measurements from table
					for(int row = 0; row < numCells; ++row){
						String rowStr = results.getRowAsString(row);
						rowStr = rowStr.substring(rowStr.indexOf(tableDelim) + 1);
						data[row + 1] += tableDelim + rowStr;
					}
					
					ResultsTable.getResultsWindow().close(false);
					WindowManager.getCurrentWindow().close();
				}
				gray.close();
				/*******************************************************************************/
				
				//Write all data to text
				String imageData = "";
				//Combine data for image
				for(int d = 0; d < data.length; ++d){
					
					String line;
					
					//Add header line on first image only
					if(d == 0 && i == 0){
						line = data[d] + "ImageName" + tableDelim + "IncImageArea" + tableDelim + "TotalImageArea" + "\n";
					}
					else{
						line = data[d] + tableDelim + color.getTitle() + tableDelim + 
								(int)(examArea + 0.5) + tableDelim + (int)(totalArea + 0.5) + "\n";
					}
					
					line = line.replaceAll(tableDelim, outDelim);
					imageData += line;
				}
				
				this.outputCellFeatures += imageData;
				
				//Pause, if appropriate
				if(this.pauseBetweenImages && i != ids.length - 1){
					IJ.showMessage("Finished with image " + color.getTitle() + ". Moving to next...");
				}
				
				//Flatten ROIs and save overlay images
				if(outDir != null && saveOverlays){
					WindowManager.setCurrentWindow(color.getWindow());
					IJ.run("Flatten");
					
					String flatPath = outDir + File.separator + "overlay-" + color.getTitle();
					ImagePlus flat = WindowManager.getCurrentImage();
					FileSaver fs = new FileSaver(flat);
					fs.saveAsTiff(flatPath);
					flat.close();
				}
				
				//Get rid of current cell ROIs
				rois.runCommand("Select All");
				rois.runCommand("Delete");
				
				IJ.showProgress((1.0*i)/ids.length);
				
				//Close window
				color.close();
			}
			
			IJ.showProgress(1.0);
			
			if(this.outDir == null || !this.writeDataToFile(this.outputCellFeatures, getCombinedOutputPath(this.outDir, "cell-features"))){
				IJ.showMessage("Error", "Couldn't write to output file. Showing in text window instead.");
				TextWindow output = new TextWindow("Cell shape & color features", this.outputCellFeatures, 750, 750);
				output.setVisible(true);
			}
			
			//Display area info for all images
			if(this.displayAreas){
				String areaString = String.format("%-30s%33s%33s\n", "Image", "Examined Area", "Total Area");
				
				for(int j = 0; j < names.length; ++j){
					areaString += String.format("%-30s%30.2f%30.2f\n", displayNames[j], examAreas[j], totalAreas[j]); 
				}
				TextWindow output = new TextWindow("Image Areas", areaString, 750, 500);
				output.setVisible(true);
			}
			
		}
		
	}
	
	/**
	 * Creates a dialog used to confirm parameter values with the user. When/if the "OK"
	 * button is clicked, the analysis begins.
	 * 
	 * @return
	 *   True if the user entered valid params and clicked OK, false otherwise
	 */
	public boolean createDialog(){
		
		boolean ready = true;
		
		GenericDialog d = new GenericDialog("Cell count setup");
		
		d.addMessage("Set the parameter values used for segmenting cells.");
		
		d.addNumericField("Size min (in pixels^2):", this.sizeMin, 2);
		d.addNumericField("Size max (in pixels^2):", this.sizeMax, 2);
		
		d.addNumericField("Circularity min:", this.circMin, 2);
		d.addNumericField("Circularity max:", this.circMax, 2);
		
		d.addMessage("Set length scale of each pixel.");
		d.addNumericField("Scale (in pixels / length unit):", this.scale, 2);
		
		d.addCheckbox("Invert current ROI selections", this.invertRoi);
		d.addCheckbox("Save cell overlay images", this.saveOverlays);
		
		d.addMessage("Color channels used to segment cells:");
		d.addCheckbox("Red", true);
		d.addCheckbox("Green", true);
		d.addCheckbox("Blue", true);
		
		d.showDialog();
		
		if(!d.wasOKed()){
			return false;
		}
		
		double params[] = new double[5];
		
		if(!getParams(d, params)){
			IJ.showMessage("Error", "Missing parameter values!");
			return false;
		}
		
		this.sizeMin = params[0];
		this.sizeMax = params[1];
		this.circMin = params[2];
		this.circMax = params[3];
		this.scale = params[4];
		
		this.invertRoi = d.getNextBoolean();
		this.saveOverlays = d.getNextBoolean();
		
		this.useChannels = new boolean[3];
		this.useChannels[0] = d.getNextBoolean();
		this.useChannels[1] = d.getNextBoolean();
		this.useChannels[2] = d.getNextBoolean();
		
		outDir = this.getOutputPath(d);
		if(outDir == null){
			ready = false;
		}
		
		return ready;
	}
	
	/**
	 * Gets the parameter values input by the user.
	 * 
	 * @param d
	 *   Dialog from which parameter values are retrieved
	 * @param params
	 *   Output destination
	 * @return
	 *   True if all values were available, false otherwise 
	 */
	public boolean getParams(GenericDialog d, double params[]){
		
		for(int i = 0; i < params.length; ++i){
			double n = d.getNextNumber();
			if(Double.isNaN(n)){
				return false;
			}
			params[i] = n;
		}
		return true;
	}
	
	/**
	 * Gets a directory chosen by the user to use as an output directory.
	 * 
	 * @return
	 *  The directory chosen by the user or null if one is not chosen.
	 */
	public String getOutputPath(GenericDialog d){
		String path = null, startPath = null;
		
		JFileChooser chooser = new JFileChooser();
		chooser.setApproveButtonText("Choose directory");
		chooser.setDialogTitle("Choose an output directory...");
		chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		
		startPath = OpenDialog.getLastDirectory();
		
		if(startPath == null){
			chooser.setCurrentDirectory(new java.io.File(OpenDialog.getDefaultDirectory()));
		}
		else{
			chooser.setCurrentDirectory(new File(startPath).getParentFile());
			chooser.setSelectedFile(new File(startPath));
		}
		
		chooser.setAcceptAllFileFilterUsed(false);
		if(chooser.showOpenDialog(d) == JFileChooser.APPROVE_OPTION){
			path = chooser.getSelectedFile().getAbsolutePath();
		}
		
		return path;
	}
	
	/**
	 * Converts the given color image to a grayscale byte image guided by the boolean values of
	 * useChannels. Before calling this method, useChannels should be set to a 3-element array with 
	 * true values indicating that the corresponding RGB channel should be included and false indicating 
	 * that it should be ignore. Included color channels are averaged to determine final pixels values.
	 * 
	 * @param color
	 *   The image to convert to grayscale
	 * @return
	 *   A grayscale, float image computed by averaging channels in the input
	 */
	public FloatProcessor convertToGrayscale(ColorProcessor color){
		
		int count = 0;
		FloatProcessor gray = new FloatProcessor(color.getWidth(), color.getHeight());
		
		for(int i = 0; i < this.useChannels.length; ++i){
			if(this.useChannels[i]){
				gray = addImages(gray, color.getChannel(i + 1, null).convertToFloatProcessor());
				++count;
			}
		}
		
		gray.multiply(1.0 / count);
		
		return gray;
	}
	
	/**
	 * Performs pixel-wise addition of the input images
	 * @param i1
	 *   First image
	 * @param i2
	 *   Second image
	 * @return
	 *   The result of i1 + i2
	 */
	public FloatProcessor addImages(FloatProcessor i1, FloatProcessor i2){
		FloatProcessor r = new FloatProcessor(i1.getWidth(), i1.getHeight());
		for(int i = 0; i < i1.getWidth() * i1.getHeight(); ++i ){
			r.setf(i, i1.getf(i) + i2.getf(i));
		}
		return r;
	}
	
	/**
	 * Performs local contrast enhancement on the current image. Relies on the CLAHE method
	 * This is built into Fiji, but requires a plugin for vanilla ImageJ. This method
	 * attempts to run this method, but fails if it can't find it in the expected places.
	 * 
	 * Find info on CLAHE here:
	 * http://rsbweb.nih.gov/ij/plugins/clahe/index.html
	 * http://fiji.sc/Enhance_Local_Contrast_(CLAHE)
	 * 
	 * @return
	 *   True if contrast adjustment was performed, false otherwise 
	 */
	public boolean adjustContrast(int blockSize){
		
		boolean ran = true;
		String cmd = null, args = null;
		
		//Determine command needed to perform CLAHE contrast adjustment
		if(this.hasCommand(this.ijContrast)){
			cmd = this.ijContrast;
			args = this.ijArgs;
		}
		else if(this.hasCommand(this.fijiContrast)){
			cmd = this.fijiContrast;
			args = this.fijiArgs;
		}
		else{
			ran = false;
		}
		
		//Run command or report it couldn't be found
		if(cmd != null){
			args = args.replace("BLOCKSIZE", "" + blockSize);
			IJ.run(cmd, args);
		}
		else{
			ran = false;
		}
		
		return ran;
	}
	
	/**
	 * Checks if command is available in the current ImageJ instance.
	 * 
	 * @param command
	 *   The command for which to search
	 * @return
	 *   True if the current ImageJ instance can recognize the command, false otherwise
	 */
	public boolean hasCommand(String command){
		
		boolean has = false;
		Hashtable table = Menus.getCommands();
		
		for(Object key : table.keySet()){
			
			if(key.toString().equals(command)){
				has = true;
				break;
			}
		}
		
		return has;
	}
	
	/**
	 * Computes Otsu threshold for the given image.
	 * 
	 * @param ip
	 *   Grayscale image for which Otsu threshold is computed. Should probably be byte-valued.
	 * @return
	 *   Threshold that maximizes interclass variance and minimizes intraclass variance
	 */
	public int getThreshold(ImageProcessor ip){

		int t, pixels;
		double pUnder, pOver, sUnder, sOver;
		double mUnder, mOver, check, max;

		int hist[] = ip.getHistogram();
		
		pixels = ip.getStatistics().pixelCount;
		pUnder = 0;
		sUnder = 0;
		sOver = 0;
		max = 0;
		t = 0;

		//Find sum of all pixel values
		for(int i = 0; i < hist.length; ++i){
			sOver += i * hist[i];
		}

		//Iterate over t values to find optimal
		for(int i = 0; i < hist.length; ++i){
			
			//Update counts of pixels under/over current t
			pUnder += hist[i];
			pOver = pixels - pUnder;
			
			//No pixels above/below t, skip this iteration
			if(pUnder == 0 || pOver == 0){
				continue;
			}
			
			//Update pixel value sums over/under current t
			sUnder += i*hist[i];
			sOver -= i*hist[i];
			
			//Get pixel means above/below current t
			mUnder = sUnder/pUnder;
			mOver = sOver/pOver;
			
			//Compute sigma and check against max
			check = (pUnder/pixels) * (pOver/pixels) * (mUnder - mOver) * (mUnder - mOver);
			if(check > max){
				max = check;
				t = i;
			}
		}
		
		return t;
	}
	
	/**
	 * Computes feature image using the given command. Applies filter (Gaussian deriv, Hessian, Laplacian)
	 * to the given image to generate feature image. The string featureCmd specifies the type and params
	 * (scale, deriv order, etc.) used in applying the filter.
	 * 
	 * Requires FeatureJ plugin (included in Fiji, url: XXX)
	 * 
	 * @param featureCmd
	 *   Command passed to FeatureJ to generate feature image 
	 * @return
	 *   Feature image corresponding to the given command/params
	 */
	public static ImageProcessor computeFeature(ImageProcessor ip, String featureCmd, String featureArgs){
		
		ImageProcessor response = null;
		
		IJ.run(featureCmd, featureArgs);
		
		return response;
	}
	
	
	/**
	 * Writes data from the given array of strings to an output text file at path. Any existing file
	 * is overwritten.
	 * 
	 * @param data
	 *   String data to write to file at path.
	 * @param path
	 *   Path of file to create.
	 * @return
	 *   True if data is successfully written, false otherwise.
	 */
	public boolean writeDataToFile(String data, String path){
		
		File outFile = new File(path);

		try{

			if(outFile.exists()){
				outFile.delete();
			}
			outFile.createNewFile();

			FileOutputStream out = new FileOutputStream(outFile);

			out.write(data.getBytes());
			
			out.flush();
			out.close();
			
		}catch (IOException e){
			return false;
		}
		
		return true;
		
	}
	
	/**
	 * Writes data from the given array of strings to an output text file at path. Any existing file
	 * is overwritten.
	 * 
	 * @param dir
	 *   Directory containing files
	 * @param files
	 *   List of files to combine
	 * @param dest
	 *   Destination to write combined file
	 * @return
	 *   True if the files are successfully combined, false otherwise
	 */
	public boolean combineDataFiles(String dir, String files[], String dest){
		
		boolean result = true;
		int count = 1;
		
		//Overwrite existing file
		File outFile = new File(dest);
		if(outFile.exists()){
			outFile.delete();
		}
		
		try {
			outFile.createNewFile();
			FileOutputStream out = new FileOutputStream(outFile);

			for(int i = 0; i < files.length; ++i){

				try{
					Scanner in = new Scanner(new File(dir + File.separator + files[i]));

					String line = in.nextLine() + "\n";

					//Keep header from first file, skip others
					if(i == 0){
						out.write(line.getBytes());
					}

					while(in.hasNext()){

						line = in.nextLine();

						line = line.substring(line.indexOf(this.outDelim));
						line = count + line + "\n";
						out.write(line.getBytes());
						++count;
					}
					
					new File(dir + File.separator + files[i]).delete();
					
				//Catch exception for this input file
				} catch (IOException e){
					result = false;
				}

			}

			out.flush();
			out.close();
		
		//Catch exception for creating/writing to output
		} catch(IOException e){
			result = false;
		}
		
		return result;
	}
	
	/**
	 * Finds unused location for outputting data. Tries to find an used name by trying to append consecutive numbers 
	 * to basename. Returns path to the first unused name.
	 *  
	 * @param dir
	 *   Directory in which output is stored.
	 * @param basename
	 *   Basename of output file, consecutive numeric suffixes ("_1", "_2", ...) are appended until an unused path
	 *   is found.
	 * @return
	 *   Unused path in dir beginning with basename
	 */
	public String getCombinedOutputPath(String dir, String basename){
		
		int counter = 1;
		boolean found = false;
		String name = null, path = null;
		
		name = basename;
		
		while(!found){
			
			path = dir + File.separator + name + ".csv";
			
			found = ! (new File(path)).exists();
			
			if(!found){
				name = basename + "_" + counter;
				counter++;
			}
		}
		
		return path;
	}
	
	/** 
	 * Attempts to create an output file and writer object for area values to be stored. Creates output
	 * file in user specified directory named "image-areas.txt" 
	 * 
	 * @return
	 *   A Formatter object to which area values can be written, null if it can't be created
	 */
	public Formatter getAreaWriter(){
		
		Formatter out = null;
		
		if(outDir != null){
			String path = outDir + File.separator + "image-areas.txt";
			
			//Try to create output file and Formatter object to write to it
			//Any errors -> return null
			try {
				
				File f = new File(path);
				if(f.exists()){
					f.delete();
				}
				
				f.createNewFile();
				
				out = new Formatter(f);
				
			} catch (IOException e) {
				out = null;
			}
		}
		
		return out;
	}
	
	/**
	 * Read default parameter values from the JAR file. The resource "parameters.config"
	 * 
	 * Never got to this...
	 * 
	 */
	public void getParameterValues(){
		Scanner in = new Scanner(RetFM_J.class.getResourceAsStream("parameters.config"));
		IJ.showMessage(in.nextLine());
	}
	
	/**
	 * Used for debugging.
	 * 
	 * @param args
	 */
	public static void main(String args[]){
		System.out.println("Shit");
	}
	
}
