#include "morphing.h"

linepairs::linepairs()
{
	max_id = 0;
}

// format of linepairs file:
// % indicates the entire line is a comment
// N (integer, indicates number of lines in the file)
// P_i P_j Q_i Q_j P_prime_i P_prime_j Q_prime_i Q_prime_j (integer VXL coordinates)
//
// Note: do not exceed 80 chars per line in the file
//

bool linepairs::load(const char* fname)
{
	vcl_ifstream infile(fname);
	char buf[80];
	int n;
	bool ok = false;

	// keep the new lines in a temporary space
	linepairs new_linepairs;

	if (!infile)
		return false;

	if (infile.good()) {
		ok = true;

		while (infile.peek() == '%')
			infile.getline(buf, 80);
		infile >> n;
		infile.getline(buf, 80);

		for (int i=0; (i<n) && (ok); i++) {
			// eliminate any comments
			while (infile.peek() == '%')
				infile.getline(buf, 80);
			// read line coordinates
			int P_i, P_j, Q_i, Q_j, P_prime_i, P_prime_j, Q_prime_i, Q_prime_j;
			if (infile.good()) {
				infile 
					>> P_i
					>> P_j
					>> Q_i
					>> Q_j
					>> P_prime_i
					>> P_prime_j
					>> Q_prime_i
					>> Q_prime_j;
				infile.getline(buf, 80);
				new_linepairs.add(P_i, P_j, Q_i, Q_j, P_prime_i, P_prime_j, Q_prime_i, Q_prime_j);
			} else 
				ok = false;
		} 
		if (ok) {
			// new line data is ok, so we copy it from the
			// temporary space to our queue, after emptying the
			// existing one
			while (pairs_.size() != 0) {
				linepair* lp;
				lp = pairs_.front();
				pairs_.pop();
				delete lp;
			}
			while (new_linepairs.pairs_.size() != 0) {
				linepair* lp;
				lp = new_linepairs.pairs_.front();
				add(lp->P_i[0], lp->P_j[0], lp->Q_i[0], lp->Q_j[0],
				    lp->P_i[1], lp->P_j[1], lp->Q_i[1], lp->Q_j[1]);
				new_linepairs.pairs_.pop();
				delete lp;
			}
			ok = true;
		} else 
			ok = false;
	} 
	return ok;
}

bool linepairs::save(const char* fname)
{
	vcl_ofstream outfile(fname);
	char buf[80];
	int n;

	if (!outfile)
		return false;

	vnl_matrix<double> P, Q, P_prime, Q_prime;
	get(P, Q, P_prime, Q_prime);

	outfile << "% Line pairs for field morphing\n";
	outfile << P.cols() << "\n";
	for (n=0; n<P.cols(); n++) 
		if (outfile.good())
			outfile << "% Line " << n << "\n"
				<< P(0,n) << " " << P(1,n) << " " 
				<< Q(0,n) << " " << Q(1,n) << " " 
				<< P_prime(0,n) << " " << P_prime(1,n) << " " 
				<< Q_prime(0,n) << " " << Q_prime(1,n) << "\n";
		else
			return false;
	
	return true;
}

void linepairs::clear()
{
	while (pairs_.size() != 0)
		pairs_.pop();
}

int linepairs::add(double P_i, double P_j, double Q_i, double Q_j, 
		            double Pp_i, double Pp_j, double Qp_i, double Qp_j)
{
	linepair* lp = new linepair;

	// we only add non-degenerate lines (ie. lines that do not have identical endpoints)
	if (((P_i == Q_i) && (P_j == Q_j)) || ((Pp_i == Qp_i) && (Pp_j == Qp_j)))
		return -1;

	lp->P_i[0] = P_i;
	lp->P_j[0] = P_j;
	lp->Q_i[0] = Q_i;
	lp->Q_j[0] = Q_j;
	lp->P_i[1] = Pp_i;
	lp->P_j[1] = Pp_j;
	lp->Q_i[1] = Qp_i;
	lp->Q_j[1] = Qp_j;
	lp->id = (max_id)++;
	pairs_.push(lp);

	return lp->id;
}

void linepairs::print(const linepair& p) {
	vcl_cerr << p.id << "\tI0: (" << p.P_i[0] << "," << p.P_j[0] << ") -> (" << p.Q_i[0] << "," << p.Q_j[0] << ")\n" 
		<< "\tI1: (" << p.P_i[1] << "," << p.P_j[0] << ") -> (" << p.Q_i[0] << "," << p.Q_j[0] << ")\n";
}

void linepairs::get(vnl_matrix<double>& P, vnl_matrix<double>& Q, 
					vnl_matrix<double>& P_prime, vnl_matrix<double>& Q_prime)
{
	linepair* lp;
	int n = pairs_.size();

	if (n > 0) {
		P.set_size(2, pairs_.size());
		Q.set_size(2, pairs_.size());
		P_prime.set_size(2, pairs_.size());
		Q_prime.set_size(2, pairs_.size());
	
		for (int i=0; i<n; i++) {
			lp = pairs_.front();
			pairs_.pop();
			P(0,i) = lp->P_i[0];
			P(1,i) = lp->P_j[0];
			Q(0,i) = lp->Q_i[0];
			Q(1,i) = lp->Q_j[0];
			P_prime(0,i) = lp->P_i[1];
			P_prime(1,i) = lp->P_j[1];
			Q_prime(0,i) = lp->Q_i[1];
			Q_prime(1,i) = lp->Q_j[1];
			pairs_.push(lp);
		}
	} else {
		P.set_size(2, 0);
		Q.set_size(2, 0);
		P_prime.set_size(2, 0);
		Q_prime.set_size(2, 0);
	}
}

bool linepairs::find_closest(int i, int j, int index, bool& isP, int& minid)
{
	linepair *lp, *bestlp;
	int n = pairs_.size();
	bool found=false;
	double dist, diP, djP, diQ, djQ, distP, distQ, mindist;

	mindist = -1;
	for (int l=0; l<n; l++) {
		lp = pairs_.front();
		pairs_.pop();
		diP = lp->P_i[index]-i;
		djP = lp->P_j[index]-j;
		diQ = lp->Q_i[index]-i;
		djQ = lp->Q_j[index]-j;
		distP = diP*diP + djP*djP;
		distQ = diQ*diQ + djQ*djQ;
		dist = vcl_min(distP, distQ);
		if ((mindist == -1) || ((mindist >=0) && 
			(dist < mindist))) {
			mindist = dist;
			minid = lp->id;
			isP = (mindist == distP);
			bestlp = lp;
		} 
		pairs_.push(lp);
	}
			
	return (mindist != -1);
}

bool linepairs::modify(int id, int index, bool isP, int ni, int nj)
{
	if ((index != 0) && (index != 1))
		return false;

	linepair* lp;
	bool found=false;
	int n = pairs_.size();

	for (int i=0; i<n; i++) {
		lp = pairs_.front();
		pairs_.pop();
		if (lp->id == id) {
			found = true;
			if (isP) {
				lp->P_i[index] = ni;
				lp->P_j[index] = nj;
			} else {
				lp->Q_i[index] = ni;
				lp->Q_j[index] = nj;
			}
		}
		pairs_.push(lp);
	}
		
	return found;
}


bool linepairs::remove(int id)
{
	bool found=false;
	int n = pairs_.size();

	for (int i=0; i<n; i++) {
		linepair* lp;

		lp = pairs_.front();
		pairs_.pop();
		if (lp->id == id) {
			found = true;
			delete lp;
		} else
			pairs_.push(lp);
	}
				
	return found;
}

linepairs linepairs::interpolate(double t)
{
	linepairs interp;
    double new_P_prime_i, new_P_prime_j, new_Q_prime_i, new_Q_prime_j;
	int n = pairs_.size();

	for (int i=0; i<n; i++) {
		linepair* lp = pairs_.front();

		pairs_.pop();
		new_P_prime_i = lp->P_i[0]*(1-t) + lp->P_i[1]*t;
		new_P_prime_j = lp->P_j[0]*(1-t) + lp->P_j[1]*t;
		new_Q_prime_i = lp->Q_i[0]*(1-t) + lp->Q_i[1]*t;
		new_Q_prime_j = lp->Q_j[0]*(1-t) + lp->Q_j[1]*t;
		pairs_.push(lp);
		interp.add(lp->P_i[0], lp->P_j[0], lp->Q_i[0], lp->Q_j[0], 
			    new_P_prime_i, new_P_prime_j, 
				new_Q_prime_i, new_Q_prime_j); 
	}

	return interp;
}

linepairs linepairs::swap()
{
	linepairs swapped;
	linepair* lp;
    double new_P_prime_i, new_P_prime_j, new_Q_prime_i, new_Q_prime_j;
	int n = pairs_.size();

	for (int i=0; i<n; i++) {
		lp = pairs_.front();
		pairs_.pop();
		swapped.add(lp->P_i[1], lp->P_j[1], lp->Q_i[1], lp->Q_j[1], 
			        lp->P_i[0], lp->P_j[0], lp->Q_i[0], lp->Q_j[0]);
		pairs_.push(lp);
	}

	return swapped;
}

linepairs linepairs::copy(int from, int to)
{
	linepairs copied;

	if ((from < 0) || (from > 1) || (to < 0) || (to > 1))
		return copied;

	linepair* lp;
    double new_P_prime_i, new_P_prime_j, new_Q_prime_i, new_Q_prime_j;
	int n = pairs_.size();

	for (int i=0; i<n; i++) {
		lp = pairs_.front();
		pairs_.pop();
		copied.add(lp->P_i[from], lp->P_j[from], lp->Q_i[from], lp->Q_j[from], 
				   lp->P_i[from], lp->P_j[from], lp->Q_i[from], lp->Q_j[from]);
		pairs_.push(lp);
	}

	return copied;
}

