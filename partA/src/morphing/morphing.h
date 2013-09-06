// DO NOT MODIFY THIS FILE EXCEPT WHERE EXCPLICTLY NOTED!!!

#ifndef _morphing_h
#define _morphing_h

#include "../vxl_includes.h"

#include "../gl/glutils.h"
#include "../imdraw/imdraw.h"

#include "linepairs.h"

// the main morphing class
class morphing {
public:
	// descriptors for all the images/input used in the algorithm
	enum im_type {I0, I1, WarpedI0, WarpedI1, Morph, Lines};
	// left/right descriptors
	enum side {Left, Right};

private:
	// the dimensions of the current image set
	int ni_;
	int nj_;
	int nplanes_;

	//
	// The images used by the morphing algorithm 
	// These images are allocated once and do not need
	// to be reallocated during morphing
	// 

	// the source image 1 (stored locally)
	vil_image_view<vil_rgb<vxl_byte> > I0_;
	// the source image 2 (stored locally)
	vil_image_view<vil_rgb<vxl_byte> > I1_;
	// the image that holds the result of warping I0
	vil_image_view<vil_rgb<vxl_byte> > warped_I0_;
	// the image that holds the result of warping I1
	vil_image_view<vil_rgb<vxl_byte> > warped_I1_;
	// the image that holds the result of morphing
	vil_image_view<vil_rgb<vxl_byte> > morph_;
	
	//
	// other major data structures used by the algorithm
	// 

	// the data structures holding the set of user-specified 
	// corresponding line pairs between image I0_ and I1_
	linepairs I0I1_linepairs_;
	// the data structures holding the set of computed corresponding
	// line pairs between image I0_ and warped_I0_
	linepairs I0W0_linepairs_;

	// the algorithm's t interpolation parameter (controls the
	// interpolation between the line pairs and the cross dissolve)
	double t_;
	// how many intermediate images to generate
	int num_images_;
	vcl_string morph_basename_;
	// should we write the warped I0 and I1 images to disc?
	bool write_warped_;
	// should we write the morphed image to disk?
	bool write_morph_;
	// the parameters of the multiple-line algorithm
	double a_;
	double b_;
	double p_;

	//
	//  Private methods
	//

	// This routine applies Beier-Neely morphing on a pair of images for
	// a given value of the t parameter and writes the results to disk
	bool morph_iteration(int iter);

	//////////////////////////////////////////////////////////////////////
	//   Methods to be implemented as Part of A.2
	//////////////////////////////////////////////////////////////////////
	
	// The main method of the feature-based metamorphosis algorithm
	// it implements the multiple-line algorithm loop on the 2nd column 
	// of page 37 of Baier & Neely's paper
	//   * source:      is the source image
	//   * linePairs:   is the pair of corresponding lines
	//   * destination: is the destination image
	//
	// This is the main function you are asked to implement in Part A.2
	// of this assignment. Your implementation should go in file
	// morphing_algorithm.cxx
	//
	static void field_warp(
			const vil_image_view<vil_rgb<vxl_byte> >& source,
			linepairs& linePairs,
			double a, double b, double p,
			vil_image_view<vil_rgb<vxl_byte> >& destination);

	
	// The second routine implements Beier-Neely morphing by calling the
	// field_warp routine
	// 
	// You are asked to implement this function as well in Part A.2
	// of the assignment. Your implementation should go in file
	// morphing_algorithm.cxx

	void compute_morph();

	//////////////////////////////////////////////////////////////////////
	//   
	//////////////////////////////////////////////////////////////////////

	// 
	// Initialization routines
	// 
	
	// routine that initializes the basic data structures of the
	// warping algorithm
	void initialize();

	// routine that sets the algorithm's default parameters
	void default_params();


	//
	// Routines & variables useful for visualizing the progress of the
	// algorithm and for algorithm debugging purposes
	//

	// is drawing enabled on the display panels?
	bool draw_enabled_;
	// should the corresponding line pairs be shown on the display panels?
	bool left_lines_visible_;
	bool right_lines_visible_;
	// which image is being displayed in each panel
	im_type left_image_;
	im_type right_image_;
	int selected_line_id_;
	bool update_panel(ImDraw *panel, bool lines_visible, im_type imt);
	// are the lines in this panel editable? returns true only if the 
	// panel is showing one of the source images
	bool editable(const ImDraw* panel);

	// PLACE ANY ADDITIONAL PRIVATE METHOD DECLARATIONS HERE 
	// AND PLACE YOUR CODE FOR THESE METHODS IN FILE 
	// morphing_algorithm.cxx

	//////////////////////////////////////////////////
	// PLACE YOUR CODE BETWEEN THESE LINES          //
	//////////////////////////////////////////////////

	//////////////////////////////////////////////////

public:
	//
    // class constructors
	// 

	// default constructor
	morphing(void);

	// constructor used when we want to show 
	// debugging information on the opengl canvas
	morphing(ImDraw* lpanel, ImDraw* rpanel);

	// 
	// methods for executing the morphing 
	// algorithm. they return false if the
	// algorithm cannot be run
	//	
	
	// run the algorithm
	bool compute();

	//
	// controlling the method's parameters
	//
	void set_a(double a);
	void set_b(double b);
	void set_p(double p);
	void set_t(double t);
	void set_num_images(int n);
	void set_morph_basename(vcl_string& str);

	// write warped images I0 and I1 to disk
	void write_warped();
	void toggle_write_warped();
	void toggle_write_morph();

	double get_a();
	double get_b();
	double get_p();
	double get_t();
	int get_num_images();
	static double get_a_default();
	static double get_b_default();
	static double get_p_default();
	static double get_t_default();
	static int get_num_images_default();

	//
	// controlling the display 
	//

	// 
	// public functions for communicating to/from the UI
	// (the associated private variables and methods are below)
	// 

	// which image to show on the specified panel 
	bool display_image(side lr, im_type imt);

	// show the lines in the specified panel
	void show_lines(side lr);
	void hide_lines(side lr);
	void toggle_lines(side lr);
	bool update_display();

	// pass pointers to the opengl display panels
	void add_panels(ImDraw* left, ImDraw* right);

	void add_line(double P_i, double P_j, double Q_i, double Q_j,  ImDraw* panel);
	void remove_line(int id);
	void modify_line(int id, bool isP, int newi, int newj, ImDraw* panel);
	void clear_lines();
	bool copy_lines(side s1, side s2);
    bool find_closest_line(int i, int j, bool& isP, int& id, const ImDraw* panel);
	int last_selected_id();

	// 
	// public functions for loading/saving linepairs
	// 
	bool save_linepairs(const char* fname);
	bool load_linepairs(const char* fname);

	// methods for setting the input to the morphing algorithm
	// the methods return true if the image im can be added to the
	// input dataset (ie. it contains image data and has the same 
	// dimensions as the already-specified images)
	bool set(im_type imt, vil_image_view<vil_rgb<vxl_byte> > im);

	// accessor functions for the inputs & outputs of the algorithm
	// the functions return an empty image if the requested information 
	// is not yet available
	bool get(im_type imt, vil_image_view<vil_rgb<vxl_byte> >& im);
	// get descriptive title of each image
	void get_title(im_type imt, vcl_string& title);

private:
	// 
	// private variables used for passing information to/from the UI
	//

	// descriptive strings for each of these images
	vcl_vector<vcl_string> im_labels;

	// flag indicating that the results have already been
	// computed
	bool morph_computed_;
	bool first_image_;

	// flag indicating that the results in the data structure
	// (if any) are out of date
	bool outdated_;

	// method that performs error checking of the user-supplied image input_im
	// and updates the variable im accordingly
	bool check_and_set_input(vil_image_view<vil_rgb<vxl_byte> > input_im,
							 vil_image_view<vil_rgb<vxl_byte> > &im);

	// we use these pointers to issue drawing commands to the
	// left and right opengl panels
	ImDraw* left_panel_;
	ImDraw* right_panel_;
};

#endif


