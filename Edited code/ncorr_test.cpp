#include "ncorr.h"

using namespace ncorr;

int disp_filewrite(DIC_analysis_input DIC_input, DIC_analysis_output DIC_output, std::string locate, std::string perspective) {
	
	//Open file to store the roi limits in each frame
	std::ofstream out_roi;
	std::string roi_file = locate + "/DataFiles/" + perspective + "/roi_limits.csv";
	out_roi.open (roi_file);
	
	for (int i = 0; i<DIC_output.disps.size(); i++) {
		Disp2D disp = DIC_output.disps[i];

		// Get the corresponding v and u Array2Ds
		const Array2D<double> &v_array = disp.get_v().get_array();
		const Array2D<double> &u_array = disp.get_u().get_array();

		// Get the corresponding ROI2D
		ROI2D disp_roi = disp.get_roi();

		//Store coordinate values of points in ROI
		std::vector<std::vector <int> > roi_pts;
		int index = 0;
		for (int p2 = 0; p2 < disp.data_width(); ++p2) {
			for (int p1 = 0; p1 < disp.data_height(); ++p1) {
				if (disp_roi(p1,p2)) {
					roi_pts.push_back(std::vector<int>());
					roi_pts[index].push_back(p1);
					roi_pts[index].push_back(p2);
					index++;
				} 
			}
		}
		
			
		//Find minimum and maximum coordinates of the ROI
		int rowmax = roi_pts[0][0], rowmin = roi_pts[0][0];
		int colmax = roi_pts[0][1], colmin = roi_pts[0][1];
	
		for (int i = 0; i < roi_pts.size(); i++) {
			//Find row min and max
			if (roi_pts[i][0] <= rowmin)
				rowmin = roi_pts[i][0];
			else if (roi_pts[i][0] >= rowmax)
				rowmax = roi_pts[i][0];
		
			//Find column min and max
			if (roi_pts[i][1] <= colmin)
				colmin = roi_pts[i][1];
			else if (roi_pts[i][1] >= colmax)
				colmax = roi_pts[i][1];
		}
		//Roi limits stored in order given below
		out_roi << rowmin << "," << rowmax << "," << colmin << "," << colmax << std::endl;
		
		//Opening files
		std::ofstream out_v;
		std::string vfile = locate + "/DataFiles/" + perspective + "/v_displacements/v_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_v.open (vfile);
		
		std::ofstream out_u;
		std::string ufile = locate + "/DataFiles/" + perspective + "/u_displacements/u_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_u.open (ufile);	
	
		//Write displacment values to file
		for (int p1 = rowmin; p1 < rowmax+1; ++p1) {
			for (int p2 = colmin; p2 < colmax+1; ++p2) {				
				if (disp_roi(p1,p2)) {
				out_v << v_array(p1,p2) << ",";
				out_u << u_array(p1,p2) << ",";
				} else {
					out_v << ",";
					out_u << ",";
				}
			} 			
			out_v << std::endl;
			out_u << std::endl;
		}
		out_v.close(); out_u.close();
		//-----------------------------------------------//
		/*
		// Trying out storing all the displacement values instead of just the roi
		//Opening files
		std::ofstream full_v;
		std::string vfull = locate + "DataFiles/" + perspective + "/v_displacements/Full_v_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_v.open (vfull);
		std::cout<<"Opened full v array" << std::endl;			
		std::ofstream full_u;
		std::string ufull = locate + "DataFiles/" + perspective + "/u_displacements/Full_u_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_u.open (ufull);
		std::cout<<"Opened full u array" << std::endl;				
		for (int p1 = 0; p1 < disp.data_height(); ++p1) {
			for (int p2 = 0; p2 < disp.data_width(); ++p2) {
				full_v << v_array(p1,p2) << ",";
				full_u << u_array(p1,p2) << ",";
			}
			full_v << std::endl;
			full_u << std::endl;
		}		
		full_v.close(); full_u.close();
		std::cout<<"Full arrays done" << std::to_string(i+1) << std::endl;		
		// End trial full storage
		//------------------------------------------------//
		*/
	}
	out_roi.close();
	/*
	save_DIC_video(locate + "/Videos/" + perspective + "/v_" + perspective + ".avi", 
				   DIC_input, 
				   DIC_output, 
				   DISP::V,
				   0.5,		// Alpha		
				   10);		// FPS

	save_DIC_video(locate + "/Videos/" + perspective + "/u_" + perspective + ".avi", 
				   DIC_input, 
				   DIC_output, 
				   DISP::U, 
				   0.5,		// Alpha
				   10);		// FPS
	*/
	return 0;
}	

/*
int strain_filewrite(strain_analysis_input strain_input, strain_analysis_output strain_output, std::string locate, std::string perspective) {
	for (int i = 0; i<strain_output.strains.size(); i++) {
		Strain2D strain = strain_output.strains[i];

		// Get the corresponding eyy, exy, and exx Array2Ds
		const Array2D<double> &eyy_array = strain.get_eyy().get_array();
		const Array2D<double> &exy_array = strain.get_exy().get_array();
		const Array2D<double> &exx_array = strain.get_exx().get_array();

		// Get the corresponding ROI2D
		ROI2D strain_roi = strain.get_roi();

		//Store coordinate values of points in ROI
		std::vector<std::vector <int> > roi_pts;
		int index = 0;
		for (int p2 = 0; p2 < strain.data_width(); ++p2) {
			for (int p1 = 0; p1 < strain.data_height(); ++p1) {
				if (strain_roi(p1,p2)) {
					roi_pts.push_back(std::vector<int>());
					roi_pts[index].push_back(p1);
					roi_pts[index].push_back(p2);
					index++;
				} 
			}
		}
	
	
		//Find minimum and maximum coordinates of the ROI
		int rowmax = roi_pts[0][0], rowmin = roi_pts[0][0];
		int colmax = roi_pts[0][1], colmin = roi_pts[0][1];
	
		for (int i = 0; i < roi_pts.size(); i++) {
			//Find row min and max
			if (roi_pts[i][0] <= rowmin)
				rowmin = roi_pts[i][0];
			else if (roi_pts[i][0] >= rowmax)
				rowmax = roi_pts[i][0];
		
			//Find column min and max
			if (roi_pts[i][1] <= colmin)
				colmin = roi_pts[i][1];
			else if (roi_pts[i][1] >= colmax)
				colmax = roi_pts[i][1];
		}
	
		//Opening files
		std::ofstream out_eyy;
		std::string eyyfile = locate + "/DataFiles/" + perspective + "/eyy_strains/eyy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_eyy.open (eyyfile);
		
		std::ofstream out_exy;
		std::string exyfile = locate + "/DataFiles/" + perspective + "/exy_strains/exy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_exy.open (exyfile);
		
		std::ofstream out_exx;
		std::string exxfile = locate + "/DataFiles/" + perspective + "/exx_strains/exx_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_exx.open (exxfile);
	
		//Write displacment values to file
		for (int p1 = rowmin; p1 < rowmax+1; ++p1) {
			for (int p2 = colmin; p2 < colmax+1; ++p2) {
				if (strain_roi(p1,p2)) {
				out_eyy << eyy_array(p1,p2) << ",";
				out_exy << exy_array(p1,p2) << ",";
				out_exx << exx_array(p1,p2) << ",";
				} else {
					out_eyy << ",";
					out_exy << ",";
					out_exx << ",";
				}
			} 
			out_eyy << std::endl;
			out_exy << std::endl;
			out_exx << std::endl;
		}

		out_eyy.close(); out_exy.close(); out_exx.close();
		
		//-----------------------------------------------//
/*		// Trying out storing all the displacement values instead of just the roi
		//Opening files
		std::ofstream full_eyy;
		std::string eyyfull = locate + "DataFiles/" + perspective + "/eyy_strains/Full_eyy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_eyy.open (eyyfull);
		std::cout<<"Opened full eyy array" << std::endl;	
		
		std::ofstream full_exy;
		std::string exyfull = locate + "DataFiles/" + perspective + "/exy_strains/Full_exy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_exy.open (exyfull);
		std::cout<<"Opened full exy array" << std::endl;	

		std::ofstream full_exx;
		std::string exxfull = locate + "DataFiles/" + perspective + "/exx_strains/Full_exx_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_exx.open (exxfull);
		std::cout<<"Opened full exy array" << std::endl;
		
		for (int p1 = 0; p1 < strain.data_height(); ++p1) {
			for (int p2 = 0; p2 < strain.data_width(); ++p2) {
				full_eyy << eyy_array(p1,p2) << ",";
				full_exy << exy_array(p1,p2) << ",";
				full_exx << exx_array(p1,p2) << ",";
			}
			full_eyy << std::endl;
			full_exy << std::endl;
			full_exx << std::endl;
		}		
		full_eyy.close(); full_exy.close(); full_exx.close();
		std::cout<<"Full strain arrays done" << std::to_string(i+1) << std::endl;		
		// End trial full storage
		//------------------------------------------------//
		*/
		/*
	}
	save_strain_video(locate + "/Videos/" + perspective + "/eyy_" + perspective + ".avi", 
					  strain_input, 
					  strain_output, 
					  STRAIN::EYY, 
					  0.5,		// Alpha
					  10);		// FPS

	save_strain_video(locate + "/Videos/" + perspective + "/exy_" + perspective + ".avi",
					  strain_input, 
					  strain_output, 
					  STRAIN::EXY, 
					  0.5,		// Alpha
					  10);		// FPS
	
	save_strain_video(locate + "/Videos/" + perspective + "/exx_" + perspective + ".avi", 
					  strain_input, 
					  strain_output, 
					  STRAIN::EXX, 
					  0.5,		// Alpha
					  10); 		// FPS
	return 0;
}*/

int dispgrad_filewrite(strain_analysis_input strain_input, strain_analysis_output strain_output, std::string locate, std::string perspective) {
	for (int i = 0; i<strain_output.strains.size(); i++) {
		Strain2D strain = strain_output.strains[i];

		// Get the corresponding eyy, exy, and exx Array2Ds
		const Array2D<double> &dvy_array = strain.get_dvy().get_array();
		const Array2D<double> &dvx_array = strain.get_dvx().get_array();
		const Array2D<double> &duy_array = strain.get_duy().get_array();
		const Array2D<double> &dux_array = strain.get_dux().get_array();

		// Get the corresponding ROI2D
		ROI2D strain_roi = strain.get_roi();

		//Store coordinate values of points in ROI
		std::vector<std::vector <int> > roi_pts;
		int index = 0;
		for (int p2 = 0; p2 < strain.data_width(); ++p2) {
			for (int p1 = 0; p1 < strain.data_height(); ++p1) {
				if (strain_roi(p1,p2)) {
					roi_pts.push_back(std::vector<int>());
					roi_pts[index].push_back(p1);
					roi_pts[index].push_back(p2);
					index++;
				} 
			}
		}
	
	
		//Find minimum and maximum coordinates of the ROI
		int rowmax = roi_pts[0][0], rowmin = roi_pts[0][0];
		int colmax = roi_pts[0][1], colmin = roi_pts[0][1];
	
		for (int i = 0; i < roi_pts.size(); i++) {
			//Find row min and max
			if (roi_pts[i][0] <= rowmin)
				rowmin = roi_pts[i][0];
			else if (roi_pts[i][0] >= rowmax)
				rowmax = roi_pts[i][0];
		
			//Find column min and max
			if (roi_pts[i][1] <= colmin)
				colmin = roi_pts[i][1];
			else if (roi_pts[i][1] >= colmax)
				colmax = roi_pts[i][1];
		}
	
		//Opening files
		std::ofstream out_dux;
		std::string duxfile = locate + "/DataFiles/" + perspective + "/dux_dispgrad/dux_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_dux.open (duxfile);
		
		std::ofstream out_duy;
		std::string duyfile = locate + "/DataFiles/" + perspective + "/duy_dispgrad/duy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_duy.open (duyfile);
		
		std::ofstream out_dvx;
		std::string dvxfile = locate + "/DataFiles/" + perspective + "/dvx_dispgrad/dvx_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_dvx.open (dvxfile);
		
		std::ofstream out_dvy;
		std::string dvyfile = locate + "/DataFiles/" + perspective + "/dvy_dispgrad/dvy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		out_dvy.open (dvyfile);
		
		
		//Write displacment values to file
		for (int p1 = rowmin; p1 < rowmax+1; ++p1) {
			for (int p2 = colmin; p2 < colmax+1; ++p2) {
				if (strain_roi(p1,p2)) {
				
				out_dux << dux_array(p1,p2) << ",";
				out_duy << duy_array(p1,p2) << ",";
				out_dvx << dvx_array(p1,p2) << ",";
				out_dvy << dvy_array(p1,p2) << ",";
				} else {
					out_dux << ",";
					out_duy << ",";
					out_dvx << ",";
					out_dvy << ",";
				}
			} 
			out_dux << std::endl;
			out_duy << std::endl;
			out_dvx << std::endl;
			out_dvy << std::endl;
		}

		out_dux.close(); out_duy.close(); out_dvx.close(); out_dvy.close();
		
		//-----------------------------------------------//
/*		// Trying out storing all the displacement values instead of just the roi
		//Opening files
		std::ofstream full_eyy;
		std::string eyyfull = locate + "DataFiles/" + perspective + "/eyy_strains/Full_eyy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_eyy.open (eyyfull);
		std::cout<<"Opened full eyy array" << std::endl;	
		
		std::ofstream full_exy;
		std::string exyfull = locate + "DataFiles/" + perspective + "/exy_strains/Full_exy_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_exy.open (exyfull);
		std::cout<<"Opened full exy array" << std::endl;	

		std::ofstream full_exx;
		std::string exxfull = locate + "DataFiles/" + perspective + "/exx_strains/Full_exx_array_" + perspective + "_frame" + std::to_string(i+1) + ".csv";
		full_exx.open (exxfull);
		std::cout<<"Opened full exy array" << std::endl;
		
		for (int p1 = 0; p1 < strain.data_height(); ++p1) {
			for (int p2 = 0; p2 < strain.data_width(); ++p2) {
				full_eyy << eyy_array(p1,p2) << ",";
				full_exy << exy_array(p1,p2) << ",";
				full_exx << exx_array(p1,p2) << ",";
			}
			full_eyy << std::endl;
			full_exy << std::endl;
			full_exx << std::endl;
		}		
		full_eyy.close(); full_exy.close(); full_exx.close();
		std::cout<<"Full strain arrays done" << std::to_string(i+1) << std::endl;		
		// End trial full storage
		//------------------------------------------------//
		*/
	}/*
	save_strain_video(locate + "/Videos/" + perspective + "/eyy_" + perspective + ".avi", 
					  strain_input, 
					  strain_output, 
					  STRAIN::EYY, 
					  0.5,		// Alpha
					  10);		// FPS

	save_strain_video(locate + "/Videos/" + perspective + "/exy_" + perspective + ".avi",
					  strain_input, 
					  strain_output, 
					  STRAIN::EXY, 
					  0.5,		// Alpha
					  10);		// FPS
	
	save_strain_video(locate + "/Videos/" + perspective + "/exx_" + perspective + ".avi", 
					  strain_input, 
					  strain_output, 
					  STRAIN::EXX, 
					  0.5,		// Alpha
					  10); 		// FPS
					  */
	return 0;
}

int main(int argc, char *argv[]) {
	/*if (argc != 2) {
		throw std::invalid_argument("Must have 1 command line input of either 'calculate' or 'dataprocess'");	
	}*/

	// Initialize DIC and strain information ---------------//
	DIC_analysis_input DIC_input;
	DIC_analysis_output DIC_output, DIC_output_Lagrangian, DIC_output_Eulerian;
	strain_analysis_input strain_input_Lagrangian, strain_input_Eulerian;
	strain_analysis_output strain_output_Lagrangian, strain_output_Eulerian;

	// Determine whether or not to perform calculations or 
	// load data (only load data if analysis has already 
	// been done and saved or else throw an exception).
	std::string input(argv[1]); //dataprocess or calculate
	std::string locate(argv[2]); // Store directory location in string variable
	std::string file_type(argv[3]); // file type for the pictures
	int file_total = std::stoi (argv[4]); // ending number for pictures
	float pixelLen = std::stof(argv[5]); // Pixels per mm in the picture
	int scale_factor = std::stoi (argv[6]); // scalefactor for the dataset
	int subregion_radius = std::stoi(argv[7]); //Subregion radius
	int threads = std::stoi(argv[8]); // # of threads
	std::string subregion_shape(argv[9]); // subregion type is circle or square
	std::string analysis_type(argv[10]); //  Decide type of analysis (small strain, large strain, discontinuous strain

	if (input == "calculate") {
		// Set images
		std::vector<Image2D> imgs;
		
		for (int i = 0; i < file_total; ++i) {
		    std::ostringstream ostr;
		    ostr << locate << "/Image_" << i << file_type;
		    imgs.push_back(ostr.str());
		}
		
		SUBREGION subregion_type;
		DIC_analysis_config DICconfig;
		// Decide shape of subregion
		if (subregion_shape == "circle") {
			subregion_type = SUBREGION::CIRCLE;
		} else if (subregion_shape == "square") {
			subregion_type = SUBREGION::SQUARE;
		}
		// Decide type of analysis (small strain, large strain, discontinuous strain
		if (analysis_type == "small") {
			DICconfig = DIC_analysis_config::NO_UPDATE;
		} else if (analysis_type == "large") {
			DICconfig = DIC_analysis_config::KEEP_MOST_POINTS;
		} else if (analysis_type == "discontinuous") {
			DICconfig = DIC_analysis_config::REMOVE_BAD_POINTS;
		}
		
		// Set DIC_input
		DIC_input = DIC_analysis_input(imgs, 							// Images
				               ROI2D(Image2D(locate + "/roi" + file_type).get_gs() > 0.5),		// ROI
					       scale_factor,                                         		// scalefactor
					       INTERP::QUINTIC_BSPLINE_PRECOMPUTE,			// Interpolation
					       subregion_type,					// Subregion shape
					       subregion_radius,                                        		// Subregion radius
					       threads,                                         		// # of threads
					       DICconfig,				// DIC configuration for reference image updates
					       true);							// Debugging enabled/disabled

		// Perform DIC_analysis    
		DIC_output = DIC_analysis(DIC_input);
		
// -----------------------------------------------------------------------------------//		
		// LAGRANGIAN
		// No process for Lagrangian (default)
		DIC_output_Lagrangian = DIC_output;

		// Set units of DIC_output_Lagrangian (provide units/pixel)
		DIC_output_Lagrangian = set_units(DIC_output_Lagrangian, "mm", pixelLen );
		
		// Set strain input (Lagrangian)
		strain_input_Lagrangian = strain_analysis_input(DIC_input,
		                                     DIC_output_Lagrangian,
		                                     subregion_type,					// Strain subregion shape
		                                     5);						// Strain subregion radius
		
		// Perform strain_analysis
		strain_output_Lagrangian = strain_analysis(strain_input_Lagrangian); 
		
// -----------------------------------------------------------------------------------//		
		// EULERIAN
		// Convert DIC_output to Eulerian perspective
		DIC_output_Eulerian = change_perspective(DIC_output, INTERP::QUINTIC_BSPLINE_PRECOMPUTE);

		// Set units of DIC_output_Eulerian (provide units/pixel)
		DIC_output_Eulerian = set_units(DIC_output_Eulerian, "mm", pixelLen);

		// Set strain input (Eulerian)
		strain_input_Eulerian = strain_analysis_input(DIC_input,
		                                     DIC_output_Eulerian,
		                                     subregion_type,					// Strain subregion shape
		                                     5);						// Strain subregion radius

		
		// Perform strain_analysis
		strain_output_Eulerian = strain_analysis(strain_input_Eulerian); 

// -----------------------------------------------------------------------------------//

		
		// Save outputs as binary
                save(DIC_input, locate + "/DataFiles/DIC_input.bin");
				
                save(DIC_output_Lagrangian, locate + "/DataFiles/DIC_output_Lagrangian.bin");
				save(DIC_output_Eulerian, locate + "/DataFiles/DIC_output_Eulerian.bin");
				
                save(strain_input_Lagrangian, locate + "/DataFiles/strain_input_Lagrangian.bin");
				save(strain_input_Eulerian, locate + "/DataFiles/strain_input_Eulerian.bin");
				
                save(strain_output_Lagrangian, locate + "/DataFiles/strain_output_Lagrangian.bin");              
                save(strain_output_Eulerian, locate + "/DataFiles/strain_output_Eulerian.bin");
				
	} else if (input == "dataprocess"){
		
		// Load inputs
		DIC_input = DIC_analysis_input::load(locate + "/DataFiles/DIC_input.bin");
		
		DIC_output_Lagrangian = DIC_analysis_output::load(locate + "/DataFiles/DIC_output_Lagrangian.bin");
		DIC_output_Eulerian = DIC_analysis_output::load(locate + "/DataFiles/DIC_output_Eulerian.bin");
		
		strain_input_Lagrangian = strain_analysis_input::load(locate + "/DataFiles/strain_input_Lagrangian.bin");
		strain_input_Eulerian = strain_analysis_input::load(locate + "/DataFiles/strain_input_Eulerian.bin");
		
		strain_output_Lagrangian = strain_analysis_output::load(locate + "/DataFiles/strain_output_Lagrangian.bin");
		strain_output_Eulerian = strain_analysis_output::load(locate + "/DataFiles/strain_output_Eulerian.bin");
		
		// Get the displacement field you want to access
		std::cout << "Number of frames is " << DIC_output_Lagrangian.disps.size() << std::endl;
		
		//DISPLACEMENTS

		disp_filewrite(DIC_input, DIC_output_Lagrangian, locate, "Lagrangian");
		disp_filewrite(DIC_input, DIC_output_Eulerian, locate, "Eulerian");

		/*//STRAINS 
		strain_filewrite(strain_input_Lagrangian, strain_output_Lagrangian, locate, "Lagrangian");
		strain_filewrite(strain_input_Eulerian, strain_output_Eulerian, locate, "Eulerian");*/
		
		//DISPLACEMENT GRADIENTS 
		dispgrad_filewrite(strain_input_Lagrangian, strain_output_Lagrangian, locate, "Lagrangian");
		dispgrad_filewrite(strain_input_Eulerian, strain_output_Eulerian, locate, "Eulerian");
		
		
	} else {
		throw std::invalid_argument("Input of " + input + " is not recognized. Must be either 'calculate' or 'dataprocess'");			
	}
	
    // Create Videos ---------------------------------------// Created in function at the start
	// Note that more inputs can be used to modify plots. 
	// If video is not saving correctly, try changing the 
	// input codec using cv::VideoWriter::fourcc(...)). Check 
	// the opencv documentation on video codecs. By default, 
	// ncorr uses cv::VideoWriter::fourcc('M','J','P','G')).

  	return 0;
}
