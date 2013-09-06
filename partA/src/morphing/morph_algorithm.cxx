// DO NOT MODIFY ANYWHERE EXCEPT WHERE EXPLICITLY NOTED!!

#include "morphing.h"

// 
// Top-level morphing routine
//
// The routine runs the basic morphing algorithm several times
// in order to generate a sequence of intermediate morphed images  
// between t=0 and t=1. 
// 
// The number of intermediate images is controlled by the 
// num_images_ parameter of the morphing class and is user-selectable
//
bool morphing::compute()
{
	int iter;
	bool ok = true;

	// if num_images > 1, we compute the t_ parameter
	// automatically before executing the morph
	if (num_images_ > 1) 
		for (iter=0; (iter<num_images_) && (ok); iter++) {
			// set the t parameter value for this iteration
			set_t((iter+1)*1.0/(num_images_+1));		
	
			// compute the morph for this setting of the parameter
			vcl_cerr << "Computing morph for t=" << get_t() << "\n";
			ok = ok && morph_iteration(iter);
		} 
	else {
		// otherwise, we just run the algorithm once, for the current
		// setting of the t parameter
		vcl_cerr << "Computing morph for t=" << get_t() << "\n";
		ok = morph_iteration(0);
	}
	return ok;
}

// 
// Top-level routine for the creation of a single morphed
// image
//
// The routine checks whether a morph can be computed
// (eg. that both images I0 and I1 have been specified),
// calls the routine that implements the Beier-Neely
// morphing algorithm, and finally writes the results to
// disk, if specified by the user
//
// the variable iter runs between 0 and num_images_ and is
// used only for file numbering when writing images to 
// disk
// 
bool morphing::morph_iteration(int iter)
{
	if ((morph_computed_ == false) || (outdated_ == true)) {

		if (outdated_ == true) {
			// we need to recalculate everything
			morph_computed_ = false;
		}

		if (((bool) I0_ == false) || 
			((bool) I1_ == false))
			// we do not have enough information yet to run the
			// algorithm, we have nothing to do
			return false;

		// ok, we have enough information to proceed

		// allocate space for all images
		warped_I0_.set_size(I0_.ni(), I0_.nj());
		warped_I1_.set_size(I0_.ni(), I0_.nj());
		morph_.set_size(I0_.ni(), I0_.nj());

		// compute the morph
		compute_morph();

		// update the compute-related flags
		morph_computed_ = true;
		outdated_ = false;

	} else {
		// the results have already been computed, so
		// we have nothing to do
	}

	// write to disk
	if (write_morph_ == true) {
		//
		// build the file name in the form <basefilename>.XXX.jpg
		// where XXX is the zero-padded iteration number
		//
		char fname[256];
		vcl_ostringstream mfname(fname);
		mfname 
			<< morph_basename_ << "." 
			<< vcl_setfill('0') << vcl_setw(3) << iter 
			<< ".jpg" << vcl_ends;
		 
		// to save, we need to access a (char *) representation
		// of the output string stream
		vcl_cerr << "writing Morph to file:" 
			<< (mfname.str()).c_str() << "\n";
		vil_save((vil_image_view<vxl_byte>)morph_, 
			     (mfname.str()).c_str());
	}
	if (write_warped_ == true) {
		char fname0[256];
		char fname1[256];
		vcl_ostringstream w0fname(fname0);
		vcl_ostringstream w1fname(fname1);

		w0fname 
			<< morph_basename_ << "." 
			<< "W0." 
			<< vcl_setfill('0') << vcl_setw(3) << iter 
			<< ".jpg" << vcl_ends;
		w1fname 
			<< morph_basename_ << "." 
			<< "W1." 
			<< vcl_setfill('0') << vcl_setw(3) << iter 
			<< ".jpg" << vcl_ends;

		vcl_cerr << "writing WarpedI0 to file " 
			<< (w0fname.str()).c_str() << "\n";
		vil_save((vil_image_view<vxl_byte>)warped_I0_, 
			     (w0fname.str()).c_str());
		vcl_cerr << "writing WarpedI1 to file " 
			<< (w1fname.str()).c_str() << "\n";
		vil_save((vil_image_view<vxl_byte>)warped_I1_, 
			     (w1fname.str()).c_str());
	}

	return true;
}

//
// Top-level implementation of the Beier-Neely morphing
// algorithm
//
// The routine should call the field_warp() routine as
// a subroutine for warping the images stored in
// variables I0_ and I1_
// 
// Specifications:
//   Inputs: The routine should assume that the following
//           private class variables contain valid data:
//
//   * The two source images, I0_ and I1_
//   * I0I1_linepairs_ holds the set of corresponding
//     line pairs between images I0_ and I1_
//   * The class variables a_, b_, p_ holding the parameters
//     of the field warping algorithm
//   * The parameter t_ that determines the in between 
//     warp between images I0_ and I1_
//
//   Outputs: The following private class variables are 
//            assumed to have valid data after the routine
//            returns:
//
//   * warped_I0_: the image holding the result of applying
//                 the field warp algorithm to image I0_
//   * warped_I1_: the image holding the result of applying
//                 the field warp algorithm to image I1_
//   * morph_: the image holding the result morphing images
//             I0_ and I1_
//   * I0W0_linepairs_: the set of corresponding line pairs
//                      between images I0_ and warped_I0_

void morphing::compute_morph()
{
	/////////////////////////////////////////
	// PLACE YOUR CODE BETWEEN THESE LINES //
	/////////////////////////////////////////

	// not so dummy any more implementation.

	linepairs lp = I0I1_linepairs_.interpolate(t_);
	linepairs lp2 = I0I1_linepairs_.swap().interpolate(1-t_);
	// just call field wrap like a baws.
	field_warp(
		I0_,
		lp,
		a_, b_, p_,
		warped_I0_);

	field_warp(
		I1_,
		lp2,
		a_, b_, p_,
		warped_I1_);
	// process for morphing.
	for (int i=0; i< ni_; i++){
		for (int j=0; j<nj_; j++) {
			morph_(i, j).r = (1-t_)*warped_I0_(i, j).r + t_*warped_I1_(i, j).r;
			morph_(i, j).g = (1-t_)*warped_I0_(i, j).g + t_*warped_I1_(i, j).g;
			morph_(i, j).b = (1-t_)*warped_I0_(i, j).b + t_*warped_I1_(i, j).b;
		}
	}


	// create a linepair where bot the original and
	// the warped lines are identical and equal to the
	// user-specified lines in image I0

	// modified linepairs.
	 I0W0_linepairs_ = lp;

	////////////////////////////////////////
}

// 
// Routine that implements the Beier-Neely field warping 
// algorithm
//
// Input variables:
//   source:      the vxl image to be warped
//   linepairs:   the set of corresponding line pairs between
//                the source image and the destination image
//   a,b,p:       the field warping parameters
//
// Output variables:
//   destination: the warped image. the routine should assume
//                that memory for this image has already been
//                allocated and that the image is of identical
//                size as the source image
// 
void morphing::field_warp(
		const vil_image_view<vil_rgb<vxl_byte> >& source, 
		linepairs& lps,
		double a, double b, double p,
		vil_image_view<vil_rgb<vxl_byte> >& destination)
{
	/////////////////////////////////////////
	// PLACE YOUR CODE BETWEEN THESE LINES //
	/////////////////////////////////////////

	/* Initialize the vectors of the P and Q lines that we will
	 * from linepars::get
	 */

	 // initialize u and v.
	 double u, v;

	 // P corresponds to the start line of the destination image. 
	 // P1 corresponds to the start line of the source image.
	 // Same applies for Q. 
	 vnl_matrix<double> P0, P1, Q0, Q1;

	 vnl_vector<double> p0, p1, q0, q1;

	 // Initialize the vectors for the source image (X1) and the dest image (X).
	 vnl_vector<double> X, X1;

	 // initialize doubles for temp result for u and v.
	 double dot_pro, mag;

	 // initialize perpendicular vectors of p and p prime.
	 vnl_vector<double> p_p, p1_p;

	 // Initialize displacement vector.
	 vnl_vector<double> disp;

	 // Initialize minimum distance.
	 double min_dist;

	 // initialize weight.
	 double w, wsum, inter_w;

	 // Initialize DSUM.
	 vnl_vector<double> dsum;

	 // initialize variables for bilinear interpolation.
	 int pix_x1, pix_x, pix_y1, pix_y, pix_final;
	 double Q11, Q12, Q21, Q22;
	 vnl_vector<double> final_pix = vnl_vector<double>(3);

	 // Initialize the dimensions of the vectors.
	 p0 = vnl_vector<double>(2);
	 p1 = vnl_vector<double>(2);
	 q0 = vnl_vector<double>(2);
	 q1 = vnl_vector<double>(2);
	 X = vnl_vector<double>(2);
	 X1 = vnl_vector<double>(2);
	 p_p = vnl_vector<double>(2);
	 p1_p = vnl_vector<double>(2);
	 disp = vnl_vector<double>(2);
	 dsum = vnl_vector<double>(2);

	 // Get the line pairs from lps.get();
	 lps.get(P1, Q1, P0, Q0);

	 for (int i=0; i<source.ni(); i++){
	 	for (int j=0; j<source.nj(); j++){

	 		// reset counters.
	 		dsum.put(0, 0);
            dsum.put(1, 0);
            wsum = 0;

	 		// initialize the destination vector.
	 		X.put(0, i);
	 		X.put(1, j);

	 		for (int l=0; l< P0.columns(); l++){

	 			// initialize p, p prime, q and q prime for each line.
	 			p0.put(0, P0.get(0, l));
	 			p0.put(1, P0.get(1, l));
	 			p1.put(0, P1.get(0, l));
	 			p1.put(1, P1.get(1, l));
	 			q0.put(0, Q0.get(0, l));
	 			q0.put(1, Q0.get(1, l));
	 			q1.put(0, Q1.get(0, l));
	 			q1.put(1, Q1.get(1, l));

	 			// Now we have sufficient information to calculate u as suggested in the paper.
	 			dot_pro = dot_product((X-p0), (q0-p0));
	 			mag = ((q0-p0).magnitude() * (q0-p0).magnitude());
	 			u = dot_pro/mag;
	 			
	 			// Now we set up the perpendicular vectors to calculate v.
	 			p_p.put(0, -(q0-p0).get(1));
	 			p_p.put(1, (q0-p0).get(0));

	 			// calculate v according to the paper.
	 			dot_pro = dot_product((X-p0), p_p);
	 			mag = (q0-p0).magnitude();
	 			v = dot_pro/mag;

	 			// Now we set up the perpendicular vector for q1 - p1.
	 			p1_p.put(0, -(q1-p1).get(1));
	 			p1_p.put(1, (q1-p1).get(0));

	 			// Now we have sufficient information to calculate X prime.
	 			X1 = p1 + u*(q1-p1) + ((v*p1_p) / (q1-p1).magnitude());

	 			// Now we find the displacement between X and X1.
	 			disp = X1 - X;

	 			// Now we need to find the shortest path as the paper suggests.
	 			if (u > 1){
	 				min_dist = (X - (q1-p1)).magnitude();
	 			}
	 			else if (u < 0) {
	 				min_dist = (X - p1).magnitude();
	 			}
	 			else {
	 				min_dist = abs(v);
	 			}
	 			inter_w = (pow((q1-p1).magnitude(), p) / (a + min_dist));
	 			w = pow(inter_w, b);
	 			dsum += (disp * w);
	 			wsum += w; 
	 		}
	 		X1 = X + (dsum/wsum);
	 		if (X1.get(0) >= destination.ni()){
	 			X1.put(0, destination.ni()-1);
	 		}
	 		if (X1.get(1) >= destination.nj()){
	 			X1.put(1, destination.nj()-1);
	 		}
	 		if (X1.get(0) < 0){
	 			X1.put(0, 0);
	 		}	 		
	 		if (X1.get(1)< 0){
	 			X1.put(1,0);
	 		}
	 		/** Perform Bilinear Interpolation **/
	 		// first of all we need to convert x and y to ints, one being of a lower value, and 
	 		// one being of a higher value. 

	 		/** Atempt at bilinear interpolation. Didn't work out for me :( I get a retarded image.

	 		pix_x = floor(X1.get(0));
	 		pix_x1 = ceil(X1.get(0));

	 		pix_y = floor(X1.get(1));
	 		pix_y1 = ceil(X1.get(1));

	 		/** Boundary Checks */
	 		/**
	 		if (pix_x1 >= destination.ni()){
	 			pix_x1 = destination.ni()-1;
	 		}
	 		if (pix_y1 >= destination.nj()){
	 			pix_y1 = destination.nj()-1;
	 		}
	 		if (pix_x < 0){
	 			pix_x = 0;
	 		}	 		
	 		if (pix_y< 0){
	 			pix_y = 0;
	 		}

	 		if (pix_x >= destination.ni()){
	 			pix_x = destination.ni()-1;
	 		}
	 		if (pix_y >= destination.nj()){
	 			pix_y = destination.nj()-1;
	 		}
	 		if (pix_x1 < 0){
	 			pix_x1 = 0;
	 		}	 		
	 		if (pix_y1< 0){
	 			pix_y1 = 0;
	 		}

	 		/* Actual Interpolation */
	 		/*
	 		Q11 = ((pix_x1-X1.get(0))*(pix_y1-X1.get(1)))/(pix_x1-pix_x*(pix_y1-pix_y));
	 		Q12 = ((X1.get(0)-pix_x1)*(pix_y1-X1.get(1)))/(pix_x1-pix_x*(pix_y1-pix_y));
	 		Q21 = ((pix_x1-X1.get(0))*(X1.get(1)-pix_y))/(pix_x1-pix_x*(pix_y1-pix_y));
	 		Q22 = ((X1.get(0)-pix_x1)*(X1.get(1)-pix_y))/(pix_x1-pix_x*(pix_y1-pix_y));

	 		final_pix.put(0, source(pix_x, pix_y).r*(Q11) + source(pix_x1, pix_y).r*(Q12) +source(pix_x, pix_y1).r*(Q21) + source(pix_x1, pix_y1).r*(Q22));
	 		final_pix.put(1, source(pix_x, pix_y).g*(Q11) + source(pix_x1, pix_y).g*(Q12) +source(pix_x, pix_y1).g*(Q21) + source(pix_x1, pix_y1).g*(Q22));
	 		final_pix.put(2, source(pix_x, pix_y).b*(Q11) + source(pix_x1, pix_y).b*(Q12) +source(pix_x, pix_y1).b*(Q21) + source(pix_x1, pix_y1).b*(Q22));

	 		destination(i,j).r = (int) final_pix.get(0);
	 		destination(i,j).g = (int) final_pix.get(1);
	 		destination(i,j).b = (int) final_pix.get(2);
		*/
	 		destination(i, j) = source(X1.get(0), X1.get(1));

	 	}
	 }
	////////////////////////////////////////
}

////////////////////////////////////////
// PLACE ANY ADDITIONAL CODE          //
// BETWEEN THESE LINES                //
////////////////////////////////////////

////////////////////////////////////////

