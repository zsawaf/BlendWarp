// DO NOT MODIFY THIS FILE!!!!

#include "morphing.h"

void morphing::add_panels(ImDraw* lpanel, ImDraw* rpanel)
{
	left_panel_ = lpanel;
	right_panel_ = rpanel;

	// pass a pointer to the panels to enable communication
	// from the panel to the morphing algorithm
	left_panel_->set_morph(this);
	right_panel_->set_morph(this);

	// since we have pointers to the display panels
	// we can enable drawing
	draw_enabled_ = true;
}

//
// Routines for interactive control of the lines used for morphing
//

// Returns true if the lines displayed over an image in one of the two
// opengl panels can be edited by the user. Only lines over images I0
// and I0 of the morphing algorithm are user-editable
//
bool morphing::editable(const ImDraw* panel)
{
	return (((panel == left_panel_) && ((left_image_ == I0) || (left_image_ == I1))) ||
		    ((panel == right_panel_) && ((right_image_ == I0) || (right_image_ == I1))));
}

// Add a pair of lines to the lineset of images I0 and I1. Initially, the lines
// over I0 and I1 are identical
void morphing::add_line(double P_i, double P_j, double Q_i, double Q_j, ImDraw* panel)
{
	// add the line to the current line set of line pairs
	if (editable(panel)) {
		outdated_ = true;

		selected_line_id_ = I0I1_linepairs_.add(P_i, P_j, Q_i, Q_j, P_i, P_j, Q_i, Q_j);
	}
	update_display();
}

// Remove the line pair of the given id from the linepair set of images I0 and I1
void morphing::remove_line(int id)
{
	outdated_ = true;
	I0I1_linepairs_.remove(id);

	update_display();
}

// Return the id of the linepair that was last edited by the user (this implements
// a limited form of user undo
int morphing::last_selected_id()
{
	return selected_line_id_;
}

// Modify the endpoint coordinates in the line pair of the given id
// See the linepairs.h file for details
void morphing::modify_line(int id, bool isP, int ni, int nj, ImDraw* panel)
{
	int index;

	if (editable(panel)) {
		if (panel == left_panel_)
			index = 0;
		else 
			index = 1;

		outdated_ = true;
		if (I0I1_linepairs_.modify(id, index, isP, ni, nj))
			selected_line_id_ = id;
	}
	update_display();
}

// Find the linepair closest to the coordinates (i,j)
// See linepairs.h for details
bool morphing::find_closest_line(int i, int j,  bool& isP, int& id, const ImDraw* panel)
{
	if (editable(panel)) {
		int index;		

		if (panel == left_panel_)
			index = 0;
		else 
			index = 1;

		return I0I1_linepairs_.find_closest(i, j, index, isP, id);
	}
	return false;
}

// Clear all the linepairs between images I0 and I1
void morphing::clear_lines()
{
	outdated_ = true;

	I0I1_linepairs_.clear();
	update_display();
}

// Copy lines from I0 to I1 or vice versa
// See linepairs.h for details
bool morphing::copy_lines(side from_side, side to_side)
{
	if ((from_side == Left) && (left_image_ == I0)) 
		if ((to_side == Right) && (right_image_ == I1))
			I0I1_linepairs_ = I0I1_linepairs_.copy(0,1);

	if ((from_side == Left) && (left_image_ == I1)) 
		if ((to_side == Right) && (right_image_ == I0))
			I0I1_linepairs_ = I0I1_linepairs_.copy(1,0);

	if ((from_side == Right) && (right_image_ == I0)) 
		if ((to_side == Left) && (left_image_ == I1))
			I0I1_linepairs_ = I0I1_linepairs_.copy(0,1);

	if ((from_side == Right) && (right_image_ == I1)) 
		if ((to_side == Left) && (left_image_ == I0))
			I0I1_linepairs_ = I0I1_linepairs_.copy(1,0);

	update_display();

	return true;
}

// Save the currently-specified line pairs between images I0 and
// I1 into a file
// See linepairs.h for details
bool morphing::save_linepairs(const char *fname)
{
	return I0I1_linepairs_.save(fname);
}

// Load line pairs from the given file 
// See linepairs.h for details
bool morphing::load_linepairs(const char *fname)
{
	if (I0I1_linepairs_.load(fname)) {
		update_display();
		outdated_ = true;
		return true;
	} else
		return false;
}

// Determine which images should be displayed in each of the
// two panels and display them, along with any lines 
bool morphing::display_image(side left_right, im_type img)
{
	bool ok = false;
	vil_image_view<vil_rgb<vxl_byte> > im;
	vcl_string title;

	if (draw_enabled_) {
		get_title(img, title);
		if (left_right == Left) {
			if (get(img, im)) 
				if (left_panel_->set(im, title)) {
					left_image_ = img;
					update_panel(left_panel_, left_lines_visible_, left_image_);
					ok = true;
				}
		} else  // we want to display an image on the right panel
			if (get(img, im)) 
				if (right_panel_->set(im, title)) {
					right_image_ = img;
					update_panel(right_panel_, right_lines_visible_, right_image_);
					ok = true;
				} 
	}
	return ok;
}

// Draw the appropriate set of lines over the specified opengl panel
bool morphing::update_panel(ImDraw *panel, bool lines_visible, im_type imt)
{
	// first erase the graphics on the panel
	panel->clear_objects();

	if (lines_visible == true) {

		switch (imt) {
		case I0:
		case I1: {
			// this is one of the input images
			// so retrieve the user-specified line pairs
			vnl_matrix<double> P0;
			vnl_matrix<double> Q0;
			vnl_matrix<double> P1;
			vnl_matrix<double> Q1;
			I0I1_linepairs_.get(P0, Q0, P1, Q1);
			// now draw the appropriate set of lines
			for (int i=0; i<P0.cols(); i++)
				if (imt == I0)
					panel->draw_object(
							new imdraw_vec(P0(1,i), P0(0,i), Q0(1,i), Q0(0,i), 
									       2.0, 3.0, 1.0, 1.0, 0.0, true));
				else 
					panel->draw_object(
							new imdraw_vec(P1(1,i), P1(0,i), Q1(1,i), Q1(0,i), 
									       2.0, 3.0, 1.0, 1.0, 0.0, true));
			}
			break;
		case WarpedI0:
		case WarpedI1:
		case Morph: {
			if (outdated_ == false) {
				// these are warped images 
				// so retrieve the interpolated line pairs
				vnl_matrix<double> P0;
				vnl_matrix<double> Q0;
				vnl_matrix<double> Pt;
				vnl_matrix<double> Qt;
				I0W0_linepairs_.get(P0, Q0, Pt, Qt);
				// now draw the appropriate set of lines
				for (int i=0; i<P0.cols(); i++)
					panel->draw_object(
							new imdraw_vec(Pt(1,i), Pt(0,i), Qt(1,i), Qt(0,i), 
										   2.0, 3.0, 1.0, 1.0, 0.0, true));
				} else
					return false;
			}
			break;
		}
	}
	return true;
}

// make the lines in the specified panel (Left or Right) visible
void morphing::show_lines(side lr)
{
	if (lr == Left)
		left_lines_visible_ = true;
	else
		right_lines_visible_ = true;

	update_display();
}

// hide the lines in the specified panel (Left or Right)
void morphing::hide_lines(side lr)
{
	if (lr == Left)
		left_lines_visible_ = false;
	else
		right_lines_visible_ = false;

	update_display();
}

// toggle visibility of lines in the specified panel
void morphing::toggle_lines(side lr)
{
	if (lr == Left)
		left_lines_visible_ = !(left_lines_visible_);
	else
		right_lines_visible_ = !(right_lines_visible_);

	update_display();
}

// refresh the line display in the two opengl panels
bool morphing::update_display()
{
	bool result = true;

	if (draw_enabled_ == true) {
		result = result && update_panel(left_panel_, left_lines_visible_, left_image_);
		result = result && update_panel(right_panel_, right_lines_visible_, right_image_);
		return result;
	} else
		return false;
}

