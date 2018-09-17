import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
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

/**
 * Computes both shape and intensity statistics of retinal cells. The cells are first extracted using a thresholding and 
 * watershed segmentation approach. This approach also includes filters that exclude cells outside of a specified size range
 * and those that are not circular enough. After segmentation, ROIs corresponding to the cells are added to the input image.
 * A set of shape and color/intensity features are computed for each cell and displayed in a results table. The set of 
 * features correspond to the complete set of measurements taken using the built-in ImageJ measurements tool. 
 * 
 * @author Mark Christopher
 *
 */
public class CellSummary_ implements PlugIn {
	
	/** Scale factor mapping pixels to length units */
	double scale = 3.15;
	
	/** Minimum size (area) of cells to count */
	double sizeMin = 5;
	
	/** Maximum size (area) of cells to count */
	double sizeMax = 125;
	
	/** Minimum circularity of cells to count */
	double circMin = 0.69;
	
	/** Maximum circularity of cells to count */
	double circMax = 1.0;
	
	/** Boolean indicating whether cells should be counted inside/outside the current ROI */
	boolean invertRoi = false;
	
	/** Boolean indicating whether the plugin should display results for each image and wait for user input */
	boolean pauseBetweenImages = false;
	
	/** Booleans indicating whether each channel should be considered */
	boolean useChannels[];
	
	/** Measurements to collect for each segmented cell */
	String measurements = "area mean standard modal min centroid center perimeter bounding fit shape " +
			"feret's integrated median skewness kurtosis area_fraction stack";
	
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
		
			//Constructor RoiManager(boolean hideWindow) never displays the window...
			RoiManager rois = new RoiManager(true);
			
			IJ.showProgress(0.0);
			
			//For every currently open image
			for(int i = 0; i < ids.length; ++i){
				
				//Get the current image/ROI
				ImagePlus p = WindowManager.getImage(ids[i]);
				
				//Remove any existing scaling from the image
				Calibration calibCopy = p.getCalibration();
				Calibration calib = new Calibration(p);
				calib.pixelWidth = 1.0;
				calib.pixelHeight = 1.0;
				p.setCalibration(calib);
				
				WindowManager.setCurrentWindow(p.getWindow());
				WindowManager.toFront(p.getWindow());
				
				String curTitle = p.getTitle();
				Roi r = p.getRoi();
				
				//Invert ROI, if needed
				if(invertRoi && r != null){
					ShapeRoi selection = new ShapeRoi(r);
					ShapeRoi all = new ShapeRoi(new Roi(0,0, p.getWidth(), p.getHeight()));
					r = all.xor(selection);
				}
				
				ImageProcessor ip = this.convertToGrayscale(p.getProcessor().convertToColorProcessor());
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
				
				//Analyze/count particles
				
				ResultsTable results = new ResultsTable();
				
				ParticleAnalyzer.setRoiManager(rois);
				int partOpts = ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES;
				ParticleAnalyzer analyze = new ParticleAnalyzer(partOpts, 0, null, 
						sizeMin*scale*scale, sizeMax*scale*scale,
						circMin, circMax);
				
				ParticleAnalyzer.setResultsTable(results);
				
				analyze.analyze(p);
				
				IJ.showMessage("Num rows + " + results.getColumn(0).length + " Rows = " + results.getRowAsString(1));
				
				IJ.run("Set Measurements...", measurements + "redirect=None decimal=3");
				
				//Take measurements of the particles
				rois.runCommand("Show All");
				rois.runCommand("Measure");
				
				if(this.pauseBetweenImages){					
					while(!IJ.spaceBarDown()){
						try{Thread.sleep(100);}
						catch(Exception ignored){}
					}
				}
				
				rois.runCommand("Select All");
				rois.runCommand("Delete");
				
				
				IJ.showProgress((1.0*i)/ids.length);
				
				//Reset scaling
				p.setCalibration(calibCopy);
			}
			
			IJ.showProgress(1.0);
		}
		
	}
	
	/**
	 * Creates a dialog used to confirm parameter values with the user. When/if the "OK"
	 * button is clicked, the analysis begans.
	 * 
	 * @return
	 *   True if the user entered valid params and clicked OK, false otherwise
	 */
	public boolean createDialog(){
		
		GenericDialog d = new GenericDialog("Cell count setup");
		
		d.addMessage("Set the parameter values used for segmenting cells.");
		
		d.addNumericField("Size min:", this.sizeMin, 2);
		d.addNumericField("Size max:", this.sizeMax, 2);
		
		d.addNumericField("Circularity min:", this.circMin, 2);
		d.addNumericField("Circularity max:", this.circMax, 2);
		
		d.addMessage("Set length scale of each pixel.");
		d.addNumericField("Scale (length units):", this.scale, 2);
		
		d.addCheckbox("Invert current ROI selections", this.invertRoi);
		d.addCheckbox("Pause between images (press space to continue)", this.pauseBetweenImages);
		
		d.addMessage("Color channels to segment cells:");
		d.addCheckbox("Red", false);
		d.addCheckbox("Green", true);
		d.addCheckbox("Blue", false);
		
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
		this.pauseBetweenImages = d.getNextBoolean();
		
		this.useChannels = new boolean[3];
		this.useChannels[0] = d.getNextBoolean();
		this.useChannels[1] = d.getNextBoolean();
		this.useChannels[2] = d.getNextBoolean();
		
		return true;
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
	 * Converts the given color image to a grayscale byte image guided by the boolean values of
	 * useChannels. Before calling this method, useChannels should be set to a 3-element array with 
	 * true values indicating that the corresponding RGB channel should be included and false indicating 
	 * that it should be ignore. Included color channels are averaged to determine final pixels values.
	 * 
	 * @param color
	 *   The image to convert to grayscale
	 * @return
	 *   A grayscale, byte image computed by averaging channels in the input
	 */
	public ImageProcessor convertToGrayscale(ColorProcessor color){
		
		int count = 0;
		FloatProcessor gray = new FloatProcessor(color.getWidth(), color.getHeight());
		
		for(int i = 0; i < this.useChannels.length; ++i){
			if(this.useChannels[i]){
				gray = addImages(gray, color.getChannel(i, null).convertToFloatProcessor());
				++count;
			}
		}
		
		gray.multiply(1.0 / count);
		
		return gray.convertToByte(false);
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
}
