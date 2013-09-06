// DO NOT MODIFY THIS FILE!!!!

#include "morphing.h"

// define a descriptive string (ie. label) for each image
// these labels are what is returned by the get_title() method

static const vcl_string emptystr("");
static const vcl_string src0str("Source I0");
static const vcl_string src1str("Source I1");
static const vcl_string warp0str("Warped I0");
static const vcl_string warp1str("Warped I1");
static const vcl_string morphstr("Morph");
static const vcl_string linesstr("Line Pairs");

// Return a descriptive label for each type of image used in the
// morphing algorithm

void morphing::get_title(im_type imt, vcl_string& str)
{
	switch (imt) {
	case I0:
		str = src0str;
		break;
	case I1:
		str = src1str;
		break;
	case WarpedI0:
		str = warp0str;
		break;
	case WarpedI1:
		str = warp1str;
		break;
	case Morph:
		str = morphstr;
		break;
	default:
	case Lines:
		str = linesstr;
	}
}

//
//  Constructor and accessor methods for setting & retrieving
//  the input/output of the morphing algorithm
//

void morphing::default_params()
{
	// initially we do not have any images loaded
	outdated_ = true;
	morph_computed_ = false;
	first_image_ = true;
	ni_ = nj_ = 0;

	// drawing is not enabled by default
	draw_enabled_ = false;
	// by default, lines are visible when drawing is enabled
	show_lines(Left); 
	show_lines(Right);
	// initially, no lines are selected
	selected_line_id_ = -1;
	// by default we don't write computed images to disk
	write_warped_ = false;
	write_morph_ = false;
	morph_basename_ = "morph";

	// set the algorithm's parameters to their default values
	a_ = get_a_default();
	b_ = get_b_default();
	p_ = get_p_default();
	// by default,compute the in-between morph between two 
	// images
	t_ = get_t_default();
	// by default, compute a single morph
	num_images_ = get_num_images_default();
}

//
// retrieving the default parameters of the morphing algorithm
//

double morphing::get_a_default()
{
	return 0.3;
}

double morphing::get_b_default()
{
	return 1.6;
}

double morphing::get_p_default()
{
	return 0.5;
}

double morphing::get_t_default()
{
	return 0.5;
}

int morphing::get_num_images_default()
{
	return 1;
}

//
// Class constructors
//

morphing::morphing() 
{
	default_params();
}

morphing::morphing(ImDraw* lpanel, ImDraw* rpanel)
{
	default_params();

	add_panels(lpanel, rpanel);
}

// 
//  Accessor methods for the class variables
// 


double morphing::get_a()
{
	return a_;
}

double morphing::get_b()
{
	return b_;
}

double morphing::get_p()
{
	return p_;
}

double morphing::get_t()
{
	return t_;
}

int morphing::get_num_images()
{
	return num_images_;
}

void morphing::set_a(double a)
{
	if (a_ != a)
		outdated_ = true;
	a_  = a;
}

void morphing::set_b(double b)
{
	if (b_ != b)
		outdated_ = true;
	b_ = b;
}

void morphing::set_p(double p)
{
	if (p_ != p)
		outdated_ = true;
	p_ = p;
}

void morphing::set_num_images(int n)
{
	if (n >= 0) 
		num_images_ = n;
}

void morphing::set_t(double t)
{
	if (t_ != t)
		outdated_ = true;
	t_ = t;
}

void morphing::set_morph_basename(vcl_string& fname)
{
	morph_basename_ = fname;
}

void morphing::write_warped()
{
	write_warped_ = true;
}

void morphing::toggle_write_warped()
{
	write_warped_ = !(write_warped_);
}

void morphing::toggle_write_morph()
{
	write_morph_ = !(write_morph_);
}

//
// Routines for setting / accessing images (mostly for communicating
// with the UI)
//

bool morphing::set(im_type imt, vil_image_view<vil_rgb<vxl_byte> > im)
{
	switch (imt) {
	case I0:
		return check_and_set_input(im, I0_);
		break;
	case I1:
		return check_and_set_input(im, I1_);
		break;
	default:
		return false;
	}
}

// methods that return an image that is the input or output
// of the morphing algorithm
// it returns an empty image if the algorithm has not been run yet on the
// input images (or if some of the images required by the algorithm
// have not been specified yet)

// this get() method is designed for returning  3-component images 
bool morphing::get(im_type imt, vil_image_view<vil_rgb<vxl_byte> >& im)
{
	switch (imt) {
	case WarpedI0:
		// we can only return a result image if the algorithm has been run
		// on the most-recently specified input images
		if ((morph_computed_ == true) &&
			(outdated_ == false)) {
			im = warped_I0_;
			return true;
		} else
			return false;
		break;
	case WarpedI1:
		// we can only return a result image if the algorithm has been run
		// on the most-recently specified input images
		if ((morph_computed_ == true) &&
			(outdated_ == false)) {
			im = warped_I1_;
			return true;
		} else
			return false;
		break;
	case Morph:
		// we can only return a result image if the algorithm has been run
		// on the most-recently specified input images
		if ((morph_computed_ == true) &&
			(outdated_ == false)) {
			im = morph_;
			return true;
		} else
			return false;
		break;
	case I0:
		if ((bool) I0_ == true) {
			im = I0_;
			return true;
		} else
			return false;
		break;
	case I1:
		if ((bool) I1_ == true) {
			im = I1_;
			return true;
		} else
			return false;
		break;
	default:
		return false;
	}
}


// routine that updates an input image of the class after first
// checking that the supplied image is valid
bool morphing::check_and_set_input(vil_image_view<vil_rgb<vxl_byte> > input_im,
									 vil_image_view<vil_rgb<vxl_byte> > &im)
{
	// does the input image contain valid image data?
	if ((bool) input_im == false)
		// if not, we return false
		return false;
		
	if (first_image_ == true) {
		// this is the first image specified, so we must set the dimensions
		// that all input images should have
		ni_ = input_im.ni();
		nj_ = input_im.nj();
		// currently we can only deal with 3-band images
		nplanes_ = 3;
		first_image_ = false;

		// add the image to the input dataset
		im = input_im;
		outdated_ = true;

		return true;
	} else
		if ((input_im.ni() == ni_) && (input_im.nj() == nj_)) {
			// ok, the image matches the dimensions of the other images
			// in the dataset, so we add it to the dataset
			im  = input_im;
			outdated_ = true;
			return true;
		} else
			// the dimensions of this image do not match the dimensions of 
			// the already-specified images so we cannot include it in the 
			// input dataset
			return false;
}



