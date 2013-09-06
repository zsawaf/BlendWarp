// DO NOT MODIFY THIS FILE EXCEPT WHERE EXPLICITLY NOTED!!!

#include "pyramid.h"

////////////////////////////////////////////////////////////////
// DO NOT MODIFY ANYTHING ABOVE THIS LINE
////////////////////////////////////////////////////////////////


bool blend(
		const vil_image_view<vxl_byte>& source0, 
		const vil_image_view<vxl_byte>& source1,
		const vil_image_view<vxl_byte>& mask,
		vil_image_view<vxl_byte>& result)
{

	// make sure image dimensions match
	if ((source0.ni() != source1.ni()) ||
		(source0.nj() != source1.nj()) ||
		(source0.nplanes() != source1.nplanes()) ||
		(source0.ni() != mask.ni()) ||
		(source0.nj() != mask.nj()))
		return false;

	// Initialize the beginning pyramids for source0, source1 and mask.
	pyramid pyr_0 = pyramid(source0);
	pyramid pyr_1 = pyramid(source1);
	pyramid pyr_M = pyramid(mask);

	// calculate the laplacians and gaussians
	vil_image_view<int> *L0 = pyr_0.L();
	vil_image_view<int> *L1 = pyr_1.L();
	vil_image_view<vxl_byte> *M = pyr_M.g();

	// calculate the total layers of the pyramid, and get the alpha value. 
	int N = pyr_0.N();
	int a = pyr_0.a();

	// allocate memory for the result matrix. 
	vil_image_view<int> S[N];

	double temp_mask;

	// Calculate the linear interprolation for each layer of the laplacians.
	for (int k=0; k<N; k++){
		// allocate memory for the kth layer of the pyramid.
		S[k] = new vil_image_view<int>(L0[k].ni(), L0[k].nj(), 3);
		for (int p=0; p<L0[k].nplanes(); p++){
			for (int rows=0; rows<L0[k].ni(); rows++){
				for (int cols=0; cols<L0[k].nj(); cols++){
					temp_mask = M[k](rows, cols)/255.0;
					S[k](rows, cols, p) = (1 - temp_mask)*L0[k](rows,cols,p) + (temp_mask*L1[k](rows,cols,p));
				}
			} 
		}
	}

	// Now we need to calculate the gaussian interprolation for the Nth layer 
	// of the pyramid.
	// set the size of the final image. 
	vil_image_view<vxl_byte> final_layer(pyr_0.g()[N].ni(), pyr_0.g()[N].nj(), 3);
	int p_limit = pyr_0.g()[N].nplanes();
	int max_rows = pyr_0.g()[N].ni();
	int max_cols = pyr_0.g()[N].nj();
	for (int p=0; p<p_limit; p++){
		for(int rows=0; rows<max_rows; rows++){
			for(int cols=0; cols<max_cols; cols++){
				temp_mask = M[N](rows,cols)/255.0;
				final_layer(rows,cols,p)= ((1 - temp_mask)*(pyr_0.g()[N](rows, cols, p))) + (temp_mask*pyr_1.g()[N](rows,cols,p));
			}
		}
	}
	// call the constructor pyramid to expand the levels we calculated in S.
	pyramid pyr_res = pyramid(S, final_layer, N, a);
	result = pyr_res.g()[0];
	return true;

}

////////////////////////////////////////////////////////////////
// DO NOT MODIFY ANYTHING BELOW THIS LINE
////////////////////////////////////////////////////////////////

bool blend2(
		const vil_image_view<vxl_byte>& source0, 
		const vil_image_view<vxl_byte>& source1,
		const vil_image_view<vxl_byte>& mask,
		vil_image_view<vxl_byte>& result)
{
	// make sure image dimensions match
	if ((source0.ni() != source1.ni()) ||
		(source0.nj() != source1.nj()) ||
		(source0.nplanes() != source1.nplanes()) ||
		(source0.ni() != mask.ni()) ||
		(source0.nj() != mask.nj()))
		return false;


	result.set_size(source0.ni(), source0.nj(), source0.nplanes());
	for (int p=0; p<source0.nplanes(); p++)
		for (int i=0; i<source0.ni(); i++)
			for (int j=0; j<source0.nj(); j++) 
				if (mask(i,j) == 0)
					result(i,j,p) = source0(i,j,p);
				else 
					result(i,j,p) = source1(i,j,p);

	return true;
}

/////////////////////////////////////////////////////
///      User interface-related routines           //
/////////////////////////////////////////////////////

// Class constructor & variable initialization
blending::blending(void) 
{

	// initially we do not have any images loaded
	outdated_ = true;
	blending_computed_ = false;
	first_image_ = true;
	draw_enabled_ = false;
	view_level_ = 0;
	view_gauss_ = true;
	view_packed_ = false;

	// initialy we display the two source images
	set_view_mode(Source0);

	source0_pyr_ = 0;
	source1_pyr_ = 0;
	blend_pyr_ = 0;
	mask_pyr_ = 0;
	N_ = 0;
}

// Top-level computation routine
// The routine just calls the blend() function after
// checking that all inputs to that function have been
// specified
bool blending::compute()
{
	if ((blending_computed_ == true) && (outdated_ == false))
		// the results have already been computed, so
		// we have nothing to do
		return true;

	if (((bool) source0_ == false) || 
		((bool) source1_ == false) || 
		((bool) mask_ == false))
		// we do not have enough information yet to run the
		// algorithm, we have nothing to do
		return false;

	// Ok, we have enough information to proceed
	
	// run the blending algorithm
	// since the algorithm operates on vxl_byte rather than
	// vil_rgb<vxl_byte> pixels, we first cast to the appropriate
	// format 
	vil_image_view<vxl_byte> result;
    blend(vil_view_as_planes(source0_), 
		  vil_view_as_planes(source1_), 
		  mask_, result);
	blend_ = result;

	// compute the pyramid of the result, for display 
	// purposes
	if (blend_pyr_)
		delete blend_pyr_;
	blend_pyr_ = new pyramid(result);


	blending_computed_ = true;
	outdated_ = false;

	return true;
}

bool blending::save_blended(const char* basename) 
{ 
     if ((blending_computed_) && (!outdated_)) { 
           vil_image_view<vxl_byte> im; 
           blend_pyr_->g(0, im); 
           char fname[256]; 
           strcpy(fname, basename); 
           strcat(fname, ".jpg"); 
           //vcl_ostringstream fname_str(fname); 
           //fname_str << basename << ".jpg"; 
           return vil_save(im, fname /*(fname_str.str()).c_str()*/); 
     } else  
           return false; 
} 

// save the results
bool blending::save_pyramid(im_type imt, bool gauss, const char* basename)
{
	char fname[256];
	vcl_ostringstream fname_str(fname);
	fname_str << basename << vcl_ends;

	pyramid* pyr;
	bool ok = false;

	switch (imt) {
	case Source0:
		if (source0_pyr_) {
			pyr = source0_pyr_;
			fname_str << ".s0" << vcl_ends;
			ok = true;
		}
		break;
	case Source1:
		if (source1_pyr_) {
			pyr = source1_pyr_;
			fname_str << ".s1" << vcl_ends;
			ok = true;
		}
		break;
	case Blend:
		if ((blending_computed_) && (!outdated_)) {
			pyr = blend_pyr_; 
			fname_str << ".b" << vcl_ends;
			ok = true;
		}
		break;
	case Mask:
		if (mask_pyr_) {
			pyr = mask_pyr_;
			fname_str << ".m" << vcl_ends;
			ok = true;
		}
		break;
	default:
		;
	}
	if (ok) 
		if (gauss) 
			pyr->dump_gauss(pyr->N(), 0, false, (fname_str.str()).c_str());
		else
			pyr->dump_laplacian(pyr->N(), 0, false, (fname_str.str()).c_str());

	return ok;
}

void blending::add_panels(ImDraw* lpanel, ImDraw* rpanel)
{
	left_panel_ = lpanel;
	right_panel_ = rpanel;

	// since we have pointers to the display panels
	// we can enable drawing
	draw_enabled_ = true;
}


// define a descriptive string (ie. label) for each image
// these labels are what is returned by the get_title() method

static const vcl_string emptystr("");
static const vcl_string src1str("Source 0");
static const vcl_string src2str("Source 1");
static const vcl_string bstr("Blended");
static const vcl_string mstr("Mask");

// Return a descriptive label for each type of image used in the
// blending algorithm

const vcl_string& blending::get_title(im_type imt)
{
	switch (imt) {
	case Source0:
		return src1str;
		break;
	case Source1:
		return src2str;
		break;
	case Blend:
		return bstr;
		break;
	default:
	case Mask:
		return mstr;
		break;
	}
}

// routine for displaying images in the left or 
// right panel
bool blending::display_images()
{
	bool ok = false;

	if (draw_enabled_) {
		ok = display_images(left_panel_, left_image_);
		ok = ok && display_images(right_panel_, right_image_);
	}
	return ok;
}

// Main display routine
//
// Controls which image is displayed in each panel
// The image to be displayed depends on the settings
// of the view_level_ and view_gauss_ variables
//    * view_level_ controls the level of the pyramid
//      that is displayed when the routine is called.
//      All levels are expanded to the level-0 image
//      before display
//    * view_gauss_ controls whether a Gauss or a
//      Laplacian pyramid level is displayed
// The function returns false if the image specified by 
// variable imt cannot be displayed
bool blending::display_images(ImDraw* panel, im_type imt)
{
	char fname[256];
	vcl_ostringstream label(fname);
	pyramid* pyr;
	bool ok = false;

	label << get_title(imt);
	if (view_packed_ == false) {
		if (view_gauss_)
			label << " Level g_l" << view_level_ << vcl_ends;
		else 
			label << " Level L_l" << view_level_ << vcl_ends;
	} else 
		if (view_gauss_)
			label << " Packed g pyramid" << vcl_ends;
		else
			label << " Packed L pyramid" << vcl_ends;

	switch (imt) {
		case Source0:
			if ((bool) source0_) {
				pyr = source0_pyr_;
				ok = true;
			}
			break;
		case Source1:
			if ((bool) source1_) {
				pyr = source1_pyr_;
				ok = true;
			}
			break;
		case Blend:
			if (blending_computed_ && (!outdated_)) {
				pyr = blend_pyr_;
				ok = true;
			}
			break;
		case Mask:
			if ((bool) mask_) {
				pyr = mask_pyr_;
				ok = true;
			}
			break;
		default:
			;
	}

	if (ok) {
		if (view_gauss_) {
			vil_image_view<vxl_byte> im2;
			if (view_packed_ == false) {
				// display a level of the Gauss pyramid
				pyr->g(view_level_, 0, im2);
				im2 = vil_crop(im2, 0, ni_, 0, nj_);
			} else 
				// display the entire pyramid
				pyr->pack_gauss(im2);
			if (im2.nplanes() == 3) {
				// convert image to RGB component format
				vil_image_view<vil_rgb<vxl_byte> > im(im2.ni(), im2.nj());
				for (int i=0; i<im2.ni(); i++)
					for (int j=0; j<im2.nj(); j++) {
						im(i,j).r = im2(i,j,0);
						im(i,j).g = im2(i,j,1);
						im(i,j).b = im2(i,j,2);
					}
				ok = panel->set(im, (label.str()).c_str());
			} else
				ok = panel->set(im2, (label.str()).c_str());
		} else {
			vil_image_view<vxl_byte> imb;
			if (view_packed_ == false) {
				// display a level of the Laplacian pyramid
				vil_image_view<int> im2;
				pyr->L(view_level_, 0, im2);
				im2 = vil_crop(im2, 0, ni_, 0, nj_);
				pyramid::int_to_ubyte(im2, imb);
			} else 
				pyr->pack_laplacian(imb);
			if (imb.nplanes() == 3) {
				vil_image_view<vil_rgb<vxl_byte> > im(imb.ni(), imb.nj());
				for (int i=0; i<imb.ni(); i++)
					for (int j=0; j<imb.nj(); j++) {
						im(i,j).r = imb(i,j,0);
						im(i,j).g = imb(i,j,1);
						im(i,j).b = imb(i,j,2);
					}
				ok = panel->set(im, (label.str()).c_str());
			} else
				ok = panel->set(imb, (label.str()).c_str());
		}

	}
	return ok;
}

// change the pyramid level being displayed in the panels
void blending::change_level(int updown)
{
	if ((N_ > 0) && (view_level_ + updown <= N_) 
		&& (view_level_ + updown >= 0)) {
		view_level_ = view_level_ + updown;
		display_images();
	}
}

// toggle between viewing Gauss vs Laplacian pyramid
void blending::toggle_view()
{
	view_gauss_ = !(view_gauss_);
	display_images();
}

// toggle between viewing packed and unpacked versions
// of the pyramids
void blending::toggle_packed()
{
	view_packed_ = (!view_packed_);
	display_images();
}

// defines which images are displayed in each viewing mode
// (there are four viewing modes)
void blending::set_view_mode(im_type imt)
{
	switch (imt) {
	case Source0: 
		view_mode_ = Source0;
		left_image_ = Source0;
		right_image_ = Source1;
		break;
	case Blend: 
		view_mode_ = Blend;
		left_image_ = Source0;
		right_image_ = Blend;
		break;
	case Source1: 
		view_mode_ = Source1;
		left_image_ = Source1;
		right_image_ = Blend;
		break;
	case Mask: 
		view_mode_ = Mask;
		left_image_ = Source0;
		right_image_ = Mask;
		break;
	}
}


// toggle between viewing source images, blended images
// or the mask
void blending::toggle_view_mode()
{
	// toggle between modes
	switch (view_mode_) {
	case Source0: 
		set_view_mode(Mask);
		break;
	case Mask: 
		set_view_mode(Blend);
		break;
	case Blend: 
		set_view_mode(Source1);
		break;
	case Source1: 
		set_view_mode(Source0);
		break;
	}
	display_images();
}

// accessor functions

bool blending::set(im_type imt, vil_image_view<vil_rgb<vxl_byte> > im)
{
	switch (imt) {
	case Source0:
		if (check_and_set_input(im, source0_)) {
			if (source0_pyr_)
				delete source0_pyr_;
			source0_pyr_ = new pyramid(vil_view_as_planes(im));
			N_ = source0_pyr_ -> N();
			set_view_mode(Source0);
			view_level_ = 0;
			view_gauss_ = true;
			display_images();
			return true;
		}
		break;
	case Source1:
		if (check_and_set_input(im, source1_)) {
			if (source1_pyr_)
				delete source1_pyr_;
			source1_pyr_ = new pyramid(vil_view_as_planes(im));
			N_ = source1_pyr_ -> N();
			set_view_mode(Source0);
			view_level_ = 0;
			view_gauss_ = true;
			display_images();
			return true;
		}
		break;
	default:
		return false;
	}
	return false;
}


bool blending::set(im_type imt, vil_image_view<vxl_byte> im)
{
	if (imt == Mask)
		if (check_and_set_input(im, mask_)) {
			if (mask_pyr_)
				delete mask_pyr_;
			mask_pyr_ = new pyramid(im);
			N_ = mask_pyr_ -> N();
			set_view_mode(Mask);
			view_level_ = 0;
			view_gauss_ = true;
			display_images();
			return true;
		}
	return false;
}

// routine that updates an input image of the class after first
// checking that the supplied image is valid
bool blending::check_and_set_input(vil_image_view<vil_rgb<vxl_byte> > input_im,
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

// routine that updates an input image of the class after first
// checking that the supplied image is valid
// identical with method above, but this version works for vxl_byte-type images
bool blending::check_and_set_input(vil_image_view<vxl_byte> input_im,
									 vil_image_view<vxl_byte> &im)
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
		first_image_ = false;

		// add the image to the input dataset
		vil_convert_cast(input_im, im);
		outdated_ = true;

		return true;
	} else
		if ((input_im.ni() == ni_) && (input_im.nj() == nj_)) {
			// ok, the image matches the dimensions of the other images
			// in the dataset, so we add it to the dataset
			im = input_im;
			outdated_ = true;
			return true;
		} else
			// the dimensions of this image do not match the dimensions of 
			// the already-specified images so we cannot include it in the 
			// input dataset
			return false;
}




