// DO NOT MODIFY THIS FILE!!!

#ifndef _pyramid_h
#define _pyramid_h

#include "../vxl_includes.h"

#include "../gl/glutils.h"
#include "../imdraw/imdraw.h"


//
// The pyramid class for implementing Laplacian &
// Gaussian pyramids
//
// The pyramid data structure stores only the Laplacian
// pyramid explicitly. The Gaussian pyramid levels are
// computed 'on the fly' as needed
//
// The Laplacian pyramid is represented as an array
// of vil_image_view<vxl_byte> images, each representing
// a different level of the pyramid. For an image of
// size (2^N + 1)*(2^N + 1), the pyramid is represented
// as an array of N-1 images, each of which is a level
// of the Laplacian pyramid, along with another vil_image_view 
// image containing the Nth level of the Gauss pyramid

class pyramid {
	// 
	// Private variables of the pyramid class
	// 

	// array of images containing the N-1 levels of the Laplacian pyramid
	vil_image_view<int>* L_;   
	// image containing the Nth level of the Gauss pyramid
	vil_image_view<vxl_byte> g_N_;
	// the total number of levels
	int N_;   
	// the "a" parameter defining the width of the 1D kernel used by
	// the expand/reduce functions
	double a_;   
	// the 1D kernel used by the expand/reduce functions
	double* w_hat_;

	// 
	// Private methods of the pyramid class
	//

	// This is the top-level routine for pyramid construction
	// It initializes the 1D kernel, builds the Gauss pyramid
	// of the image passed as a parameter, and finally builds
	// the Laplacian pyramid. 
	void build(const vil_image_view<vxl_byte>& im);

	// 
	// The reduce() and expand() functions. You will
	// have to implement both these functions
	//
	// The input in both functions is an image of dimension 
	// (2^k + 1)x(2^k + 1) with 1<=k<=N as well as the 1D
	// kernel w_hat to be used for smoothing/interpolation
	//
	static void reduce(const vil_image_view<vxl_byte> im,
		               const double* w_hat,
		               vil_image_view<vxl_byte>& im_red);
	static void expand(const vil_image_view<vxl_byte> im, 
		               const double* w_hat,
					   vil_image_view<vxl_byte>& im_exp);
	// this is identical to the expand() routine above, but it operates
	// on images with pixels of type int
	static void expand(const vil_image_view<int> im, 
		               const double* w_hat,
					   vil_image_view<int>& im_exp);

    // Pyramid-packing functions. These functions take an array of images as
	// input, each of which is 1/2 the size of the previous one, and 'packs'
	// them into a single image for visualization purposes
	void pack(vil_image_view<vxl_byte>* pyr, int l, int i, int j, 
              vil_image_view<vxl_byte>& im) const;
	// Same routine but implemented for int-pixel images (used for packing
	// Laplacian levels, which can have negative values and thereforecannot be
	// represented as vxl_byte images
	void pack(vil_image_view<int>* pyr, int l, int i, int j, 
              vil_image_view<int>& im) const;
public:
	//
	// pyramid constructors
	//

	// constructing the pyramid from a supplied image
	pyramid(const vil_image_view<vxl_byte>& im);
	// constructor with the kernel's a-parameter specified explicitly
	pyramid(const vil_image_view<vxl_byte>& im, double a);

	// basic accessor functions
	int N() const;
	double a() const;

	//
	// Construct a pyramid data structure from an array of N-1 Laplacian 
	// levels and the N-th level of the Gauss pyramid.
	//  
	// You will likely use this consructor in your implementation of the
	// blend() function
	// 
	pyramid(const vil_image_view<int>* L, const vil_image_view<vxl_byte> g_N, int N, double a);

	// 
	// Methods for accessing levels of the Gauss or the Laplacian pyramid. 
	// The Gauss pyramid versions of these functions compute the requested 
	// level on the fly
	// 

	// Store level l of the Gauss pyramid in image g_l
	// The routine returns false if l is outside the valid range
	bool g(int l, vil_image_view<vxl_byte>& g_l) const; 

	// Store level l of the Laplacian pyramid in image L_l
	// The routine returns fale if l is outside the valid range
	bool L(int l, vil_image_view<int>& L_l) const;

	// Compute level l1 of the Gaussian pyramid and then 
	// expand it to the size of level l2. 
	// The routine returns false if l1,l2 are outside the
	// valid range, or if l2 > l1.
	bool g(int l1, int l2, vil_image_view<vxl_byte>& g_l) const; 

	// Compute level l1 of the Laplacian pyramid and then 
	// expand it to the size of level l2. 
	// The routine returns false if l1,l2 are outside the
	// valid range, or if l2 > l1.
	bool L(int l1, int l2, vil_image_view<int>& L_l) const; 

	// Return an array that holds the entire Gauss pyramid
	// The routine reconstructs the Gauss pyramid from the Laplacian
	// pyramid images
	vil_image_view<vxl_byte>* g() const;

	// Return an array that holds the entire Laplacian pyramid
	vil_image_view<int>* L() const;
	                     
	//
	// Utility functions 
	//

	// Store in variable im an image that contains a 'packed' 
	// version of the Gauss pyramid
	void pack_gauss(vil_image_view<vxl_byte>& im) const;

	// Store in variable im an image that contains a 'packed' 
	// version of the Laplacian pyramid
	void pack_laplacian(vil_image_view<vxl_byte>& im) const;

	// convert an int-type image (used for storing signed Laplacian levels)
    // to a ubyte image
	static void int_to_ubyte(const vil_image_view<int>& imi, vil_image_view<vxl_byte>& imb);

	// Crop and pad an image so that it becomes square and has 
	// size (2^N+1)x(2^N+1),
	void crop_to_power_of_2plus1(
		    const vil_image_view<vxl_byte>& im, 
		    vil_image_view<vxl_byte>& im_crop); 

	// 
	// Write each level of the Laplacian pyramid as a separate 
	// image on disk. If expand_to_l2=true, the levels are expanded
	// so that the image size is equal to that of images at level l2
	//
	// The routine also outputs an image containing a packed version
	// of the entire pyramid
	// 
	// Filenames are of the form <basename>.L.x.jpg where x is the level
	// or                        <basename>.L.x_exp.jpg 
	// when the expand_to_l2 flag is true
	// and                       <basename>.L.pack.jpg for the packed image
	//
	// Use these routines for debugging your pyramid construction
	// code, as they allow you to access individual levels of the
	// pyramid
	void dump_laplacian(int l1, int l2, bool expand_to_l2, const char* basename);
    // The analogous method for dumping gaussian pyramid levels
	void dump_gauss(int l1, int l2, bool expand_to_l2, const char* basename);

};


// 
// Top-level routine that implements pyramid blending
//
// This routine should be written by you
// 
// Input:
//     source0, source1: the two images to be blended
//                       they must have identical size
//                       but do not have be of size 
//                       (2^N+1)x(2^N+1)
//     blending mask:    pixel values are typically either
//                       0 or 255 (ie the mask is usually binary)
// Output:
//     result:           the result of the pyramid
//                       blending operation
//     return value:     returns true if blending can be performed
//                       (eg. all input images are of the same size, etc)
//
static bool blend(const vil_image_view<vxl_byte>& source0, 
		             const vil_image_view<vxl_byte>& source1,
		             const vil_image_view<vxl_byte>& mask,
		             vil_image_view<vxl_byte>& result);

// Dummy blending routine that constructs the result
// by copying pixels from source0(i,j) if mask(i,j)=0
// and copies them from source1(i,j) otherwise
static bool blend2(const vil_image_view<vxl_byte>& source0, 
		              const vil_image_view<vxl_byte>& source1,
		              const vil_image_view<vxl_byte>& mask,
		              vil_image_view<vxl_byte>& result);



//////////////////////////////////////////////////////////////
//                 The blending class                       //
//////////////////////////////////////////////////////////////

//
// The blending class provides the interface between the 
// blend() routine above and the viscomp user interface
// 
// Your implementation should not depend on or use ANYTHING 
// in this class, so you can completely ignore the definitions
// below

class blending {
public:
	// descriptors for all the images/input used in the algorithm
	enum im_type {Source0, Source1, Blend, Mask};

private:
	// the dimensions of the current image set
	int ni_;
	int nj_;
	int N_;
	//int nplanes_;
	bool first_image_;

	// pointers to pyramid data structures for each of the 
	// locally-stored images
	pyramid *source0_pyr_;
	pyramid *source1_pyr_;
	pyramid *blend_pyr_;
	pyramid *mask_pyr_;

	// display control flags
	int view_level_;
	bool view_gauss_;
	bool view_packed_;
	im_type view_mode_;

	im_type left_image_;
	im_type right_image_;
    bool draw_enabled_;

	ImDraw* left_panel_;
	ImDraw* right_panel_;

	// the first source image (stored locally)
	vil_image_view<vil_rgb<vxl_byte> > source0_;
	// the second source image (stored locally)
	vil_image_view<vil_rgb<vxl_byte> > source1_;
	// the mask image (stored locally)
	vil_image_view<vxl_byte> mask_;
	// the blended image 
	vil_image_view<vil_rgb<vxl_byte> > blend_;

public:
	// routines for saving image pyramids
	//   imt      specifies the type of pyramid to be saved
	//   gauss    specifies whether to output the gauss or the laplacian pyramid
	//   basename the base filename of the result
	// the routine returns false if the save operation failed
	bool save_pyramid(im_type imt, bool gauss, const char* basename);
	// save the result of the blending operation
	// the routine returns false if the save operation failed
	bool save_blended(const char* basename);

	
	// default constructor
	blending(void);

	// run the blending algorithm on a set of
	// previously-specified input images
	bool compute();

	// pass pointers to the opengl display panels
	void add_panels(ImDraw* left, ImDraw* right);

	// change the level of the pyramid displayed in the
	// opengl panel
	void change_level(int updown);
	// toggle the display of Gaussian or Laplacian levels
	void toggle_view();
	// toggle between viewing individual levels or a packed
	// version of the entire pyramid
	void toggle_packed();
	// toggle the display of source images, blending result
	// and mask
	void toggle_view_mode();
	// specify one of the four viewing modes of the interface:
	//   Source0:  displays Source0 on left and Source1 on right
	//   Mask:              Source0 on left and Mask    on right
	//   Blend:             Source0             Blend 
	//   Source1:           Source1             Blend
	void set_view_mode(im_type imt);


	// 
	// public functions for communicating to/from the UI
	// (the associated private variables and methods are below)
	// 

	// methods for setting the input to the blending algorithm
	// the methods return true if the image im can be added to the
	// input dataset (ie. it contains image data and has the same 
	// dimensions as the already-specified images)
	bool set(im_type imt, vil_image_view<vil_rgb<vxl_byte> > im);
	bool set(im_type imt, vil_image_view<vxl_byte> im);

	// get descriptive title of each image
	const vcl_string& get_title(im_type imt);

private:

	bool display_images();
	bool display_images(ImDraw* panel, im_type imt);

	// descriptive strings for each of these images
	vcl_vector<vcl_string> im_labels_;

	// flag indicating that the results have already been
	// computed
	bool blending_computed_;
	// flag indicating that the results in the data structure
	// (if any) are out of date
	bool outdated_;

	// method that performs error checking of the user-supplied image input_im
	// and updates the variable im accordingly
	bool check_and_set_input(vil_image_view<vil_rgb<vxl_byte> > input_im,
							 vil_image_view<vil_rgb<vxl_byte> > &im);
	bool check_and_set_input(vil_image_view<vxl_byte> input_im,
							 vil_image_view<vxl_byte> &im);
};


#endif


