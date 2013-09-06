// DO NOT MODIFY THIS FILE 

#ifndef _linepairs_h
#define _linepairs_h

#include "../vxl_includes.h"

///////////////////////////////////////////////////////
// A class & methods for manipulating pairs of lines //
///////////////////////////////////////////////////////


//
// Structure defining a pair of corresponding line segments
// in two images (say, I0 and I1)
//
//  { (Pi[0],Pj[0])->(Qi[0],Qj[0]) , (Pi[1],Pj[1])->(Qi[1],Qj[1]) }
// 
// where the first line segment is in the first image and the second
// line segment is in the second image, (Pi[k],Pj[k]) is the line 
// segment's first endpoint and (Qi[k],Qj[k]) is the line segment's
// second endpoint
//

typedef struct linepair_struct {
	// Coordinates of the endpoints of the line segment
	double P_i[2], P_j[2], Q_i[2], Q_j[2];
	// The line pair's id. Each pair of corresponding line segments
	// is given a unique id that is used to identify it for 
	// interactive editing/deletion operations
	int id;
} linepair;

// 
// The linepairs class
//
// This class models the set of user-specified corresponding
// line pairs in a pair of images (usually I0 and I1).
// 
// The class is implemented as a queue of pointers to 
// linepair structures
//
class linepairs {
	// the maximum id of lines in the line set; this is used
	// for generating new unique line pair id's
	int max_id;
	vcl_queue<linepair *> pairs_;
public:
	linepairs();
	// Add a pair of lines to the existing set and return the id 
	// assigned to that pair
	int add(double P_i, double P_j, double Q_i, double Q_j,
		     double Pp_i, double Pp_j, double Qp_i, double Qp_j);
	// Modify the coordinates of an endpoint in a line pair:
	//  id:    the id of the linepair 
	//  index: 0 if the endpoint is for the line on image I0 and
	//         1 if the endpoint is for the line on image I1
	//  isP:   true if it is the P endpoint and false if it is the Q endpoint
	//  ni,nj: the new coordinates of the endpoint
	// The routine returns true if the line pair was found in the set and
	// modified and false if the line pair with the given id was not found
	bool modify(int id, int index, bool isP, int ni, int nj);
	// remove the line pair of given id from the set of line pairs
	// the routine returns true if the line pair was found and removed
	bool remove(int id);
	// find the linepair whose endpoint is closest to image coordinates (i,j)
	//  index: 0 if the coordinates refer to image I0 and 1 if they refer to I1
	//  id:    the id of the line pair whose endpoint is closest to (i,j)
	//  isP:   true if that endpoint is the P endpoint and if it is the Q endpoint
	bool find_closest(int i, int j, int index, bool& isP, int& id);

	//
	// linepair interpolation routine
	//
	// NOTE: you will probably want to use this routine in your implementation
	//       of the compute_morph() method
	//
	// Given an interpolation 
	// parameter t, the routine returns a new set of linepairs such that
	// for every pair of corresponding lines in the current line pair set
	//      { (Pi[0],Pj[0])->(Qi[0],Qj[0]) , (Pi[1],Pj[1])->(Qi[1],Qj[1]) }
	// the new set contains the pair
	//      { (Pi[0],Pj[0])->(Qi[0],Qj[0]) , (Ai,Aj)->(Bi,Bj) }
	// where
	//      Ai = (1-t)*Pi[0] + t*Pi[1], Aj = (1-t)*Pj[0] + t*Pj[1]
	//      Bi = (1-t)*Qi[0] + t*Qi[1], Bj = (1-t)*Qj[0] + t*Qj[1]
	//
	linepairs interpolate(double t);

	//
	// linepair swapping routine
	//
	// NOTE: you will probably want to use this routine in your implementation
	//       of the compute_morph() method
	//
	// The routine returns a new line pair set where the order of the line
	// segments has been reversed. Specifically, for every line pair
	//  { (Ai,Aj)->(Bi,Bj) , (Ci,Cj)->(Di,Dj) }
	// the new set will contain the pair
	//  { (Ci,Cj)->(Di,Dj) , (Ai,Aj)->(Bi,Bj) }
	linepairs swap();

	// 
	// routine that copies lines from I0 to I1 and vice versa
	// 
	// The routine returns a new line pair set where each line pair
	// 
	// { (Pi[0],Pj[0])->(Qi[0],Qj[0]) , (Pi[1],Pj[1])->(Qi[1],Qj[1]) }
	// 
	// is transformed to the pair
	//
	// { (Pi[from],Pj[from])->(Qi[from],Qj[from]) , (Pi[from],Pj[from])->(Qi[from],Qj[from]) }
    // 
	// where the input parameter from is either 0 or 1
	//
	linepairs copy(int from, int to);

	// dumping linepair data into a matrix
	//
	// NOTE: you will probably want to use this routine in your implementation
	//       of the field_warp() method
	//
	// copy all the line pair coordinates into four 2xL matrices where
	// L is the number of lines in the line set and the matrices are as
	// follows:
	// 
    //   P0 holds the (i,j) coordinates of the P endpoints of lines in image I0
	//   Q0                                    Q                             I0
	//   P1                                    P                             I1
	//   Q1                                    Q                             I1
	void get(vnl_matrix<double>& P0, vnl_matrix<double>& Q0, 
		     vnl_matrix<double>& P1, vnl_matrix<double>& Q1);

	// delete all line pairs, producing an empty line pair set
	void clear();

	// Load linepair data from a file of the given name
	// A linepair file has the following structure:
	// 
	//
    // L        (integer, indicates number of line pairs in the file)
    // P0_i P0_j Q0_i Q0_j P1_i P1_j Q1_i Q1_j (integer image coordinates)
	// 
	// % any line starting with a '%' character is treated as a comment
	//
    //
    // Note: do not exceed 80 chars per line in the file
    //
	bool load(const char* fname);

	// save linepair data into a file of the given name
	bool save(const char* fname);

	// print the linepair data to stderr
	static void print(const linepair& lp);
};

#endif

