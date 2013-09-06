// DO NOT MODIFY THIS FILE EXCEPT WHERE EXPLICITLY NOTED!!!

#include "pyramid.h"

////////////////////////////////////////////////////////////////
//          The pyramid class constructor routines            //
////////////////////////////////////////////////////////////////


pyramid::pyramid(const vil_image_view<vxl_byte>& im)
{
	// set the default value for the alpha parameter to 0.4
	a_ = 0.4;

	build(im);
}

pyramid::pyramid(const vil_image_view<vxl_byte>& im, double a)
{
	a_ = a;

	build(im);
}

//
// This is the main routine for building the Gauss and the
// Laplacian pyramid
//
// It takes as input an image (of arbitrary dimensions MxP)
// and builds a pyramid whose 0 level is of size (2^N+1)x(2^N+1)
// and N is the smallest power of 2 such that 2^N+1 > max(M,P)
//
void pyramid::build(const vil_image_view<vxl_byte>& im)
{
	int ni, nj;
	int l;
	int n;
	int i,j, plane;

	//
	// initialize the smoothing kernel
    //

	// the kernel always has length 5 pixels
	w_hat_ = new double[5];
	// shift pointer by two so that kernel indices are in
	// the range [-2,2]
	w_hat_ += 2;
	// set the values of the kernel
	w_hat_[2] = w_hat_[-2] = 0.25 - a_/2;
	w_hat_[1] = w_hat_[-1] = 0.25;
	w_hat_[0] = a_;

	// 
	// Allocate the arrays holding the Gauss and Laplacian pyramid images
	// 

	// compute the number of levels, N_
	ni = (int)ceil(log((double)im.ni()-1)/log(2));
	nj = (int)ceil(log((double)im.nj()-1)/log(2));
	N_ = vcl_max(ni,nj);

	// 
    // a pyramid is represented as an array of N_+1 images, 0,...,N_
	// where the n-th image is of dimension (2^n+1)x(2^n+1)
	// 

	// allocate a temporary array holding the Gaussian pyramid images
	vil_image_view<vxl_byte>* g_temp = new vil_image_view<vxl_byte>[N_+1];

	// allocate the array holding the Laplacian pyramid images
	L_ = new vil_image_view<int>[N_];

	// 
	// Build the pyramids
	//

	// Step 1: Copy the original image into level 0 of the Gauss pyramid

	// allocate space and initialize the level-0 image with zeros, since 
	// the input image may not occupy all (2^N_+1)x(2^N+1) pixels
	g_temp[0].set_size((int)pow(2,N_)+1, (int)pow(2,N_)+1, im.nplanes());
	g_temp[0].fill(0);

	// copy the image so that pixel im(0,0) goes to pixel g_temp[0](0,0)
	vil_copy_to_window(im, g_temp[0], 0, 0);


	// Step 2: Compute levels 1,...,N_ of the Gauss pyramid

	for (l=1; l<=N_; l++) 
		// each level is a reduced version of the image immediately below it
		reduce(g_temp[l-1], w_hat_, g_temp[l]);

	// Step 3: Compute levels 0,...,N_-1 of the Laplacian pyramid

	for (l=0; l<N_; l++) {
		// allocate space for the l-th level of the Laplacian pyramid
		L_[l] = vil_image_view<int>(g_temp[l].ni(), g_temp[l].nj(), 
			                              g_temp[l].nplanes());
		vil_image_view<vxl_byte> gtmp;

		// each level is computed by the difference
		//   L_[l] = g_[l] - expand(g_[l+1])
		// i.e., it represents all the "details" in g_[l] that are
		// "lost" when the image is reduced from g_[l] to g_[l+1]

		expand(g_temp[l+1], w_hat_, gtmp);
		vil_math_image_difference(g_temp[l], gtmp, L_[l]);
	}

	// Store just the N_ level image of the Gauss pyramid
	g_N_ = g_temp[N_];
}

// Create a pyramid data structure from a sequence of
// Laplacian pyramid levels and the N-th level Gauss image
pyramid::pyramid(const vil_image_view<int>* L, 
				 vil_image_view<vxl_byte> g_N, 
				 int N, double a)
{
	int l;

	N_ = N;
	a_ = a;

	// the kernel always has length 5 pixels
	w_hat_ = new double[5];
	// shift pointer by two so that kernel indices are in
	// the range [-2,2]
	w_hat_ += 2;
	// set the values of the kernel
	w_hat_[2] = w_hat_[-2] = 0.25 - a_/2;
	w_hat_[1] = w_hat_[-1] = 0.25;
	w_hat_[0] = a_;

	L_ = new vil_image_view<int>[N_];
	for (l=0; l<N; l++)
		vil_copy_deep(L[l], L_[l]);

	vil_copy_deep(g_N, g_N_);
}

//
// Pyramid accessor routines
// 


// Return the total number of levels in the pyramid
int pyramid::N() const
{
	return N_;
}

// Return the alpha value of the w_hat kernel
double pyramid::a() const
{
	return a_;
}

// Copy level l of the Laplacian pyramid into the supplied
// image; the routine returns FALSE if (l is not in the 
// range [0,...,N_-1]
bool pyramid::L(int l, vil_image_view<int>& L_l) const
{
	if ((l >=0) && (l < N_)) {
		vil_copy_deep(L_[l], L_l);
		return true;
	} else
		return false;
}

// Copy level l1 of the Laplacian pyramid into the supplied
// image, and expend it to level l2; the routine returns 
// FALSE if (l1 is not in the  range [0,...,N_-1]) or
// if (l1 < l2)
bool pyramid::L(int l1, int l2, vil_image_view<int>& L_l) const
{
	if ((l1 >=0) && (l1 < N_) && (l1 >= l2)) {
		vil_image_view<int> temp1, temp2;

		// get level l1 of the pyramid
		L(l1, temp1);

		// expand the level to level l2
		for (int l=l1; l>l2; l--) {
			expand(temp1, w_hat_, temp2);
			temp1 = temp2;
		}
		vil_copy_deep(temp1, L_l);

		return true;
	} else
		return false;
}

// Copy all levels of the Laplacian pyramid into an array of
// images and return a pointer to that array

vil_image_view<int>* pyramid::L() const
{
	// allocate the array holding the images
	vil_image_view<int>* L_temp = new vil_image_view<int>[N_];

	// copy the corresponding levels into the array
	for (int l=0; l<N_; l++) 
		vil_copy_deep(L_[l], L_temp[l]);
	
	return L_temp;
}


////////////////////////////////////////////////////////////////
//        Routines for computing the Gauss pyramid            //
//                                                            //
//   The gauss pyramid is not stored explicitly, so we        //
//   need to reconstruct it from the Laplacian images upon    //
//   request                                                  //
//                                                            //
////////////////////////////////////////////////////////////////

// a helper function that adds two images and also performs intensity bound
// checking to avoid overflows when assigning numbers to unsigned byte
// images
static vil_image_view<vxl_byte> add_g_and_L(vil_image_view<vxl_byte>& g, 
						                    vil_image_view<int>& L)
{
	vil_image_view<vxl_byte> result(g.ni(), g.nj(), g.nplanes());
	
	for (int p=0; p<g.nplanes(); p++)
		for (int i=0; i<g.ni(); i++)
			for (int j=0; j<g.nj(); j++) {
				int value = g(i,j,p) + L(i,j,p);
				result(i,j,p) = vcl_min(vcl_max(value, 0), 255);
			}

	return result;
}


// Return an array that holds the entire Gauss pyramid
// The routine reconstructs the Gauss pyramid from the Laplacian
// pyramid images

vil_image_view<vxl_byte>* pyramid::g() const
{
	// allocate the array holding the images
	vil_image_view<vxl_byte>* g_temp = new vil_image_view<vxl_byte>[N_+1];

	// the N_-th level of the Gauss pyramid is only level stored in the
	// pyramid data structure
	g_temp[N_] = g_N_;

	// the remaining levels have to be computed
	for (int l=N_; l>0; l--) {
		vil_image_view<vxl_byte> gtmp;
		// compute the l-th level of the Gauss pyramid using the
		// formula g_(l-1) = expand(g_l) + L_[l-1]
		expand(g_temp[l], w_hat_, gtmp);
		g_temp[l-1] = add_g_and_L(gtmp, L_[l-1]);
	}

	return g_temp;
}


// Compute the l-th level of the Gauss pyramid from the 
// Laplacian pyramid images; the routine
// returns FALSE if l is not in the range [0,...,N_]
bool pyramid::g(int l, vil_image_view<vxl_byte>& g_l) const
{
	return g(l, l, g_l);
}

// Compute the l1-th level of the Gauss pyramid from the 
// Laplacian pyramid images. The routine then expands 
// the l1-th level image repeatedly so that its size 
// becomes (2^l2 + 1)x(2^l2 + 1).
//
// The routine returns FALSE if l1,l2 are outside the 
// range [0,...,N_] and/or if l2 > l1.
// 
// If (l1=l2), the routine simply computes the l1-level image of the
// Gauss pyramid


bool pyramid::g(int l1, int l2, vil_image_view<vxl_byte>& g_l) const
{
	int l;
	vil_image_view<vxl_byte> temp;

    if ((l1 >= 0) && (l1 <= N_) &&
		(l2 >= 0) && (l2 <= N_) &&
		(l2 <= l1)) {

        // the N_-th level of the Gauss pyramid is only level stored in the
        // pyramid data structure
		vil_image_view<vxl_byte> g_l1;
		g_l1 = g_N_;

		// the remaining levels have to be computed
		for (l=N_; l>l1; l--) {
			vil_image_view<vxl_byte> gtmp;
			// compute the l-th level of the Gauss pyramid using the
			// formula g_(l-1) = expand(g_l) + L_[l-1]
			expand(g_l1, w_hat_, gtmp);
			temp = add_g_and_L(gtmp, L_[l-1]);
			g_l1 = temp;

		}

		// create gl_2 by expanding until level l2
		for (l=l1, g_l=g_l1; l>l2; l--) {
			// implement the formula g_l = expand(g_l)			
			expand(g_l, w_hat_, temp);
			g_l = temp;
		}
		return true;
	} else
		return false;
}

////////////////////////////////////////////////////////////////
// DO NOT MODIFY ANYTHING ABOVE THIS LINE
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//               The expand/reduce routines                   //
////////////////////////////////////////////////////////////////

//
// The REDUCE() routine
// 

// Given image of dimensions (2^N+1)x(2^N+1), return an image
// of dimensions (2^(N-1)+1)x(2^(N-1)+1) that is a smoothed and
// subsampled version of the original
// 
// Parameters:
//    im:     the input image
//    w_hat:  the 5-pixel smoothing kernel (w_hat[-2],...,w_hat[2])
//
// Boundary treatment: When the kernel extends beyond the image 
//                     boundary, you should compute the convolution
//                     (1) by taking into account only the kernel
//                     elements that overlap with the image
//                     and (2) by dividing the result by the sum of 
//                     these kernel elements

void pyramid::reduce(const vil_image_view<vxl_byte> im, 
			         const double* w_hat,
			         vil_image_view<vxl_byte>& im_red)
{

	
	double normal_first = 1 / (w_hat[0] + w_hat[1] + w_hat[2]);
	double normal_last = 1 / (w_hat[-2] + w_hat[-1] + w_hat[0]);	
	double red_normal = 1 / (w_hat[-2] + w_hat[-1] + w_hat[0] + w_hat[1] + w_hat[2]);
	// set the size of the reduced image to half of the original image.
	im_red.set_size((im.ni()/2)+1, (im.nj()/2)+1, im.nplanes());
	// create a temp image, and set the size of its width to half.
	vil_image_view<vxl_byte> red_width(im.ni(), (im.nj()/2)+1, im.nplanes());

	// counter for keeping track where to store the pixel in the reduced images.
	int pix_place;

	double result;

	for(int p=0; p < im.nplanes(); p++){
		for(int rows=0; rows < im.ni(); rows++){
			// reset the pixel place counter.
			pix_place = 0;
			for(int cols=0; cols<im.nj(); cols+=2){
				/* Handle Boundary Cases */
				// First boundary case is when cols is on the left edge of the image. 
				// then we have to multiply the positive w_hat with the corresponding pixels.
				if (cols == 0) {
					result = normal_first *
					((w_hat[0]*im(rows, 0, p)) +
					 (w_hat[1]*im(rows, 1, p)) +
					 (w_hat[2]*im(rows, 2, p)));
					red_width(rows, pix_place, p) = (int) result;
				}
				else if (cols == im.nj() - 1){
					red_width(rows, pix_place, p) =
					normal_last *
					((w_hat[-2]*im(rows, cols - 2, p)) +
					 (w_hat[-1]*im(rows, cols - 1, p)) +
					 (w_hat[0]*im(rows, cols, p)));
					red_width(rows, pix_place, p) = (int) result;
				}
				else {
					result = red_normal *
					((w_hat[-2]*im(rows, cols - 2, p)) +
					 (w_hat[-1]*im(rows, cols - 1, p)) +
					 (w_hat[0]*im(rows, cols, p)) +
					 (w_hat[1]*im(rows, cols + 1, p)) +
					 (w_hat[2]*im(rows, cols + 2, p)));
					red_width(rows, pix_place, p) = (int) result;

				} 
				pix_place ++;
			}
		}
	}
	pix_place = 0;
	for (int p=0; p<im.nplanes(); p++){
		for(int cols=0; cols<red_width.nj(); cols++){
			// reset pixel counter.
			pix_place=0;
			for(int rows=0; rows<red_width.ni(); rows+=2){
				/* same as previous cases */
				if(rows == 0){
					result = normal_first *
					((w_hat[0]*red_width(0, cols, p)) +
					 (w_hat[1]*red_width(1, cols, p)) +
					 (w_hat[2]*red_width(2, cols, p)));
					im_red(pix_place, cols, p) = (int) result;
				}
				else if (rows == red_width.ni() - 1){
					result = normal_last *
					((w_hat[-2]*red_width(rows-2, cols, p)) +
					 (w_hat[-1]*red_width(rows-1, cols, p)) +
					 (w_hat[0]*red_width(rows, cols, p)));
					im_red(pix_place, cols, p) = (int) result;

				}
				else {
					result = red_normal *
					((w_hat[-2]*red_width(rows-2, cols, p)) +
					 (w_hat[-1]*red_width(rows-1, cols, p)) +
					 (w_hat[0]*red_width(rows, cols, p)) +
					 (w_hat[1]*red_width(rows+1, cols, p)) +
					 (w_hat[2]*red_width(rows+2, cols, p)));
					im_red(pix_place, cols, p) = (int) result;

				}
				pix_place++;
			}
		}
	}
}

//
// The EXPAND() routine
// 

// Given image of dimensions (2^(N-1)+1)x(2^(N-1)+1), return an image
// of dimensions (2^N+1)x(2^N+1) that is an interpolated version
// of the original using the w_hat as the interpolating kernel
// 
// Parameters:
//    im:     the input image
//    w_hat:  the 5-pixel smoothing kernel (w_hat[-2],...,w_hat[2])
//
//
// Boundary treatment: When the kernel extends beyond the image 
//                     boundary, you should compute the convolution
//                     (1) by taking into account only the kernel
//                     elements that overlap with the image
//                     and (2) by dividing the result by the sum of 
//                     these kernel elements
void pyramid::expand(const vil_image_view<vxl_byte> im, 
		             const double* w_hat, 
		             vil_image_view<vxl_byte>& im_exp)
{

	// set the size for the expanded image.
	im_exp.set_size((2*im.ni())-1, (2*im.nj())-1, im.nplanes());
	// create a new temp image expanded by columns.
	vil_image_view<vxl_byte> temp_expand(im.ni(), (2*im.nj())-1, im.nplanes());
	// temp value to store intermediate results
	double result;
	// normalized w_hat.
	/* Keep in mind that I purposely reweighted the images the way they are now */
	double normalize = 1 / (w_hat[-2] + w_hat[0]);
	// normalize the expand
	double odd_normal = 1 / (w_hat[-1] + w_hat[1]);
	double even_normal = 1 / (w_hat[-2] + w_hat[0] + w_hat[2]);
	for(int p=0; p < im_exp.nplanes(); p++){
		for(int rows=0; rows<im_exp.ni()/2; rows++){
			for(int cols=0; cols<im_exp.nj(); cols++){
				/* Handle boundary cases on both edges */
				// In this case, we need to get the center pixel
				// and the 2nd pixel, and multiply them by the corresponding
				// kernel.
				if (cols == 0){
					result = normalize *
					 (w_hat[0] * im(rows, 0, p) +
					 (w_hat[2] * im(rows, 1, p)));
					temp_expand(rows, 0, p) = (int) result;
				}
				else if (cols == im_exp.nj()-1) {
					result = normalize *
					(w_hat[-2] * im(rows, im.nj()-2, p) +
					(w_hat[0] * im(rows, im.nj()-1, p)));
					temp_expand(rows, cols, p) = (int) result;
				}
				// now we handle the middle cases, there are two different cases,
				// one case is for even pixels, and one case is for odd pixels. 

				// In the case of even pixels, we are going to take the center pixel,
				// the center-2 pixel, and the center + 2 pixel and multiply them by
				// the corresponding weights.
				else if (cols%2 == 0) {
					result = even_normal *
					((w_hat[-2] * im(rows, (cols/2)-1, p) +
					(w_hat[0] * im(rows, cols/2, p)) +
					(w_hat[2] * im(rows, (cols/2)+1, p))));
					temp_expand(rows, cols, p) = (int) result;
				}
				// In the case of odd pixels, the center pixel does not exist. Therefore,
				// we take the center-1, and the center+1 and multiply them by the 
				// corresponding weights.
				else {
					result = odd_normal *
					((w_hat[-1] * im(rows, cols/2, p)) +
					 (w_hat[1] * im(rows, (cols/2)+1, p)));
					temp_expand(rows, cols, p) = (int) result;
				}
			}
		}
	}
	// right now we need to expand along the columns. 
	for(int p=0; p < im_exp.nplanes(); p++){
		for(int cols=0; cols<im_exp.nj(); cols++){
			for(int rows=0; rows<im_exp.ni(); rows++){
				/* Handle boundary cases on both edges */
				// In this case, we need to get the center pixel
				// and the 2nd pixel, and multiply them by the corresponding
				// kernel.
				if (rows == 0){
					result = normalize *
					 (w_hat[0] * temp_expand(0, cols, p) +
					 (w_hat[2] * temp_expand(1, cols, p)));
					im_exp(0, cols, p) = (int) result;
				}
				else if (rows == im_exp.nj()-1) {
					result = normalize *
					(w_hat[-2] * temp_expand(temp_expand.ni()-2, cols, p) +
					(w_hat[0] * temp_expand(temp_expand.ni()-1, cols, p)));
					im_exp(rows, cols, p) = (int) result;
				}
				// now we handle the middle cases, there are two different cases,
				// one case is for even pixels, and one case is for odd pixels. 

				// In the case of even pixels, we are going to take the center pixel,
				// the center-2 pixel, and the center + 2 pixel and multiply them by
				// the corresponding weights.
				else if (rows%2 == 0) {
					result = even_normal *
					((w_hat[-2] * temp_expand((rows/2)-1, cols, p) +
					(w_hat[0] * temp_expand(rows/2, cols, p)) +
					(w_hat[2] * temp_expand((rows/2)+1, cols, p))));
					im_exp(rows, cols, p) = (int) result;
				}
				// In the case of odd pixels, the center pixel does not exist. Therefore,
				// we take the center-1, and the center+1 and multiply them by the 
				// corresponding weights.
				else {
					result = odd_normal *
					((w_hat[-1] * temp_expand(rows/2, cols, p)) +
					 (w_hat[1] * temp_expand((rows/2)+1, cols, p)));
					im_exp(rows, cols, p) = (int) result;
				}
			}
		}
	}
}

void pyramid::expand(const vil_image_view<int> im, 
		             const double* w_hat, 
		             vil_image_view<int>& im_exp)
{
		// set the size for the expanded image.
	im_exp.set_size((2*im.ni())-1, (2*im.nj())-1, im.nplanes());
	// create a new temp image expanded by columns.
	vil_image_view<int> temp_expand(im.ni(), (2*im.nj())-1, im.nplanes());
	// temp value to store intermediate results
	double result;
	// normalized w_hat.
	/* Keep in mind that I purposely reweighted the images the way they are now */
	double normalize = 1 / (w_hat[-2] + w_hat[0]);
	// normalize the expand
	double odd_normal = 1 / (w_hat[-1] + w_hat[1]);
	double even_normal = 1 / (w_hat[-2] + w_hat[0] + w_hat[2]);
	for(int p=0; p < im_exp.nplanes(); p++){
		for(int rows=0; rows<im_exp.ni()/2; rows++){
			for(int cols=0; cols<im_exp.nj(); cols++){
				/* Handle boundary cases on both edges */
				// In this case, we need to get the center pixel
				// and the 2nd pixel, and multiply them by the corresponding
				// kernel.
				if (cols == 0){
					result = normalize *
					 (w_hat[0] * im(rows, 0, p) +
					 (w_hat[2] * im(rows, 1, p)));
					temp_expand(rows, 0, p) = (int) result;
				}
				else if (cols == im_exp.nj()-1) {
					result = normalize *
					(w_hat[-2] * im(rows, im.nj()-2, p) +
					(w_hat[0] * im(rows, im.nj()-1, p)));
					temp_expand(rows, cols, p) = (int) result;
				}
				// now we handle the middle cases, there are two different cases,
				// one case is for even pixels, and one case is for odd pixels. 

				// In the case of even pixels, we are going to take the center pixel,
				// the center-2 pixel, and the center + 2 pixel and multiply them by
				// the corresponding weights.
				else if (cols%2 == 0) {
					result = even_normal *
					((w_hat[-2] * im(rows, (cols/2)-1, p) +
					(w_hat[0] * im(rows, cols/2, p)) +
					(w_hat[2] * im(rows, (cols/2)+1, p))));
					temp_expand(rows, cols, p) = (int) result;
				}
				// In the case of odd pixels, the center pixel does not exist. Therefore,
				// we take the center-1, and the center+1 and multiply them by the 
				// corresponding weights.
				else {
					result = odd_normal *
					((w_hat[-1] * im(rows, cols/2, p)) +
					 (w_hat[1] * im(rows, (cols/2)+1, p)));
					temp_expand(rows, cols, p) = (int) result;
				}
			}
		}
	}
	// right now we need to expand along the columns. 
	for(int p=0; p < im_exp.nplanes(); p++){
		for(int cols=0; cols<im_exp.nj(); cols++){
			for(int rows=0; rows<im_exp.ni(); rows++){
				/* Handle boundary cases on both edges */
				// In this case, we need to get the center pixel
				// and the 2nd pixel, and multiply them by the corresponding
				// kernel.
				if (rows == 0){
					result = normalize *
					 (w_hat[0] * temp_expand(0, cols, p) +
					 (w_hat[2] * temp_expand(1, cols, p)));
					im_exp(0, cols, p) = (int) result;
				}
				else if (rows == im_exp.nj()-1) {
					result = normalize *
					(w_hat[-2] * temp_expand(temp_expand.ni()-2, cols, p) +
					(w_hat[0] * temp_expand(temp_expand.ni()-1, cols, p)));
					im_exp(rows, cols, p) = (int) result;
				}
				// now we handle the middle cases, there are two different cases,
				// one case is for even pixels, and one case is for odd pixels. 

				// In the case of even pixels, we are going to take the center pixel,
				// the center-2 pixel, and the center + 2 pixel and multiply them by
				// the corresponding weights.
				else if (rows%2 == 0) {
					result = even_normal *
					((w_hat[-2] * temp_expand((rows/2)-1, cols, p) +
					(w_hat[0] * temp_expand(rows/2, cols, p)) +
					(w_hat[2] * temp_expand((rows/2)+1, cols, p))));
					im_exp(rows, cols, p) = (int) result;
				}
				// In the case of odd pixels, the center pixel does not exist. Therefore,
				// we take the center-1, and the center+1 and multiply them by the 
				// corresponding weights.
				else {
					result = odd_normal *
					((w_hat[-1] * temp_expand(rows/2, cols, p)) +
					 (w_hat[1] * temp_expand((rows/2)+1, cols, p)));
					im_exp(rows, cols, p) = (int) result;
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////
// DO NOT MODIFY ANYTHING BELOW THIS LINE
////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////
//           Pyramid class utility routines                   //
////////////////////////////////////////////////////////////////

//
// Crop an input image so that it becomes square,
// has dimensions 2^N + 1, and 2^N+1 is both smaller
// and closest to the largest dimension of the image
//

void pyramid::crop_to_power_of_2plus1(
		const vil_image_view<vxl_byte>& im,
		vil_image_view<vxl_byte>& im_crop)
{
	int ni, nj, N;

	// compute the number of levels, N
	ni = (int)floor(log((int)im.ni())/log(2));
	nj = (int)floor(log((int)im.nj())/log(2));
	N = vcl_max(ni,nj);

	im_crop = vil_crop(im, 0, (int)(pow(2,N)+1), 0, (int)(pow(2,N)+1));
}

// 
// Write levels l1 to l2 of the Gaussian pyramid to 
// disk. If expand_to_l2=true, all images are expanded 
// to the size of the level l2 images. 
//
// The routine also outputs a packed version of the pyramid
// into a single image
void pyramid::dump_gauss(int l1, int l2, bool expand_to_l2, const char* basename)
{
	int i;

	for (int l=l2; l<=l1; l++) {
		vil_image_view<vxl_byte> im;
		char fname[256];
		vcl_ostringstream pyr_fname(fname);

		if (expand_to_l2 == true) {
			pyr_fname << basename << ".g_exp." << l << ".jpg" << vcl_ends;
			g(l, l2, im);
		} else {
			pyr_fname << basename << ".g." << l << ".jpg" << vcl_ends;
			g(l, im);
		}

		 
		// to save, we need to access a (char *) representation
		// of the output string stream
		vcl_cerr 
			<< "Saving Gauss pyramid level " << l << " to file "
			<< (pyr_fname.str()).c_str() << "\n";
		vil_save(im, (pyr_fname.str()).c_str());
	}

	// save a packed version of the pyramid
	char fname[256];
	vcl_ostringstream pack_fname(fname);
	vil_image_view<vxl_byte> im;

	pack_fname << basename << ".g.pack.jpg" << vcl_ends;
	vcl_cerr 
		<< "Saving packed Gauss pyramid to file "
		<< (pack_fname.str()).c_str() << "\n";
	pack_gauss(im);
	vil_save(im, (pack_fname.str()).c_str());
}

void pyramid::dump_laplacian(int l1, int l2, bool expand_to_l2, const char* basename)
{
	int i;

	if (l1 >= N_)
		l1 = N_-1;
	if (l1 < 0) 
		l1 = 0;
	if (l2 >= N_)
		l2 = N_;
	if (l2 < 0)
		l2 = 0;

	for (int l=l2; l<=l1; l++) {
		vil_image_view<int> im;
		vil_image_view<vxl_byte> imb;
		char fname[256];
		vcl_ostringstream pyr_fname(fname);

		if (expand_to_l2 == true) {
			pyr_fname << basename << ".L_exp." << l << ".jpg" << vcl_ends;
			L(l, l2, im);
		} else {
			pyr_fname << basename << ".L." << l << ".jpg" << vcl_ends;
			L(l, im);
		}

		 
		// to save, we need to access a (char *) representation
		// of the output string stream
		vcl_cerr 
			<< "Saving Laplacian pyramid level " << l << " to file "
			<< (pyr_fname.str()).c_str() << "\n";
		int_to_ubyte(im, imb);
		vil_save(imb, (pyr_fname.str()).c_str());
	}

	// save a packed version of the pyramid
	char fname[256];
	vcl_ostringstream pack_fname(fname);
	vil_image_view<vxl_byte> im;

	pack_fname << basename << ".L.pack.jpg" << vcl_ends;
	vcl_cerr 
		<< "Saving packed Laplacian pyramid to file "
		<< (pack_fname.str()).c_str() << "\n";
	pack_laplacian(im);
	vil_save(im, (pack_fname.str()).c_str());

}

// 
// Return an image that contains all levels of the Laplacian pyramid
// This routine is useful for pyramid visualization purposes
// 
void pyramid::pack_laplacian(vil_image_view<vxl_byte>& imb) const
{
	// allocate space for the output image; the packed image has
	// 1/2 more columns and the same number of rows as the original image 
	vil_image_view<int> im(2*(L_[0].ni()-1), 
		                   L_[0].nj(), 
						   L_[0].nplanes());
	// Note: a more efficient packing is possible:
	// allocate space for the output image; the packed image has
	// 1/2 more columns and the same number of rows as the original image 
	//vil_image_view<int> im(L_[0].ni()+(L_[0].ni()-1)/2+1, 
	//	                   L_[0].nj(), 
	//					   L_[0].nplanes());

	// fill it with zeros
	im.fill(0);

	// call the recursive packing routine
	pack(L_, N_, 0, 0, im);

	// since the Laplacian images contain signed byte values, we first convert it to an
	// unsigned byte image by mapping the range [-127..128] to [0..255]
	int_to_ubyte(im, imb);
}

// convert an int-type image (used for storing signed Laplacian images)
// to a ubyte image
void pyramid::int_to_ubyte(const vil_image_view<int>& imi, vil_image_view<vxl_byte>& imb)
{
	imb.set_size(imi.ni(), imi.nj(), imi.nplanes());
	for (int p=0; p<imi.nplanes(); p++)
		for (int i=0; i<imi.ni(); i++)
			for (int j=0; j<imi.nj(); j++)
				imb(i,j,p) = vcl_min<int>(imi(i,j,p)+127,255);

}

// 
// Return an image that contains all levels of the Gauss pyramid
// This routine is useful for pyramid visualization purposes. 
//
void pyramid::pack_gauss(vil_image_view<vxl_byte>& im) const
{
	// Since the Gaussian pyramid is not stored explicitly, we first
	// compute & store the Gaussian pyramid in a temporary array
	// of images
	vil_image_view<vxl_byte>* g_temp = g();

	// allocate space for the output image; the packed image has
	// the same number of rows as the original image and twice the 
	// columns
	im.set_size(2*(g_temp[0].ni()-1), g_temp[0].nj(), g_temp[0].nplanes());
	// note: a tigher packing is also possible
	//im.set_size(g_temp[0].ni()+(g_temp[0].ni()-1)/2+1, 
	//	                        g_temp[0].nj(), 
	//                           g_temp[0].nplanes());

	// fill it with zeros
	im.fill(0);
	// call the recursive pyramid-packing routine
	pack(g_temp, N_, 0, 0, im);
}

// Recursive packing routine for unsigned byte images
void pyramid::pack(vil_image_view<vxl_byte>* pyr, int l, int i, int j, 
				   vil_image_view<vxl_byte>& im) const
{
	if (l > 0) {
		vil_copy_to_window(pyr[0], im, i, j);
		if (l > 1) {
			vil_copy_to_window(pyr[1], im, i+pyr[0].ni(), j);
			if (l > 2) {
				vil_copy_to_window(pyr[2], im, i+pyr[0].ni(), j+pyr[1].nj());
				if (l >= 3) 
					pack(pyr+3, l-3, i+pyr[0].ni()+pyr[2].ni(), j+pyr[1].nj(), im);
			}
		}
	}
}

// Recursive packing routine for signed byte images (useful for 
// packing Laplacian images). 
void pyramid::pack(vil_image_view<int>* pyr, int l, int i, int j, 
				   vil_image_view<int>& im) const
{
	if (l > 0) {
		vil_copy_to_window(pyr[0], im, i, j);
		if (l > 1) {
			vil_copy_to_window(pyr[1], im, i+pyr[0].ni(), j);
			if (l > 2) {
				vil_copy_to_window(pyr[2], im, i+pyr[0].ni(), j+pyr[1].nj());
				if (l >= 3) 
					pack(pyr+3, l-3, i+pyr[0].ni()+pyr[2].ni(), j+pyr[1].nj(), im);
			}
		}
	}
}


