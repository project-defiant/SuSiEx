/*
 * data.hpp
 *
 *  Created on: Dec 27, 2021
 *      Author: kyuan
 */

#ifndef DATA_HPP_
#define DATA_HPP_

# include <iostream>
# include <vector>
# include <cmath>
# include <cstdlib>
# include <cstring>

# define LDTYPE float

extern int npop;

class softpar
{
public:
	softpar() : out_dir(""), out_name(""), chr(""), start(0), end(0), level(0.95), \
	min_purity(0.5), pth(1e-5), tol(1e-4), mult_step(false), n_sig(5), max_iter(100), nthreads(1), \
	key_by(2) {}

	std::vector<int> n_gwas;
	std::string out_dir, out_name;
	std::string chr;
	long start, end;
	double level, min_purity, pth, tol;
	bool mult_step;
	int n_sig, max_iter, nthreads;

	/*
	 * 0 rsID, 1 position ref / alt, 2 position + rsID
	 * 2 as default
	 */
	int key_by;
};

class varCoord
{
public:
	std::string id, a1, a2;
	long pos;
};

inline bool cmp_rs(const varCoord& a, const varCoord& b)
{
	return a.id < b.id;
}

inline bool eql_rs(const varCoord& a, const varCoord& b)
{
	return a.id == b.id;
}

inline bool cmp_pos_ref(const varCoord& a, const varCoord& b)
{
	if(a.pos != b.pos)
		return a.pos < b.pos;
	if(a.a1 != b.a1)
		return a.a1 < b.a1;
	return a.a2 < b.a2;
}

inline bool eql_pos_ref(const varCoord& a, const varCoord& b)
{
	return a.pos == b.pos && a.a1 == b.a1 && a.a2 == b.a2;
}

inline bool cmp_pos_rs(const varCoord& a, const varCoord& b)
{
	if(a.pos != b.pos)
		return a.pos < b.pos;
	return a.id < b.id;
}

inline bool eql_pos_rs(const varCoord& a, const varCoord& b)
{
	return a.pos == b.pos && a.id == b.id;
}

static bool (*var_cmp_func)(const varCoord& a, const varCoord& b) = cmp_pos_rs;

static bool (*var_eql_func)(const varCoord& a, const varCoord& b) = eql_pos_rs;

class sumstats
{
public:
	varCoord coord;
	int idx;
	double beta;
	double pval, log10p;
};

inline bool operator<(const sumstats& a, const sumstats& b)
{
	return var_cmp_func(a.coord, b.coord);
}

inline bool operator==(const sumstats& a, const sumstats& b)
{
	return var_eql_func(a.coord, b.coord);
}

class ldref
{
public:
	varCoord coord;
	int idx;
	double frq;
};

inline bool operator<(const ldref& a, const ldref& b)
{
	return var_cmp_func(a.coord, b.coord);
}

inline bool operator==(const ldref& a, const ldref& b)
{
	return var_eql_func(a.coord, b.coord);
}

inline int stat_amb(const varCoord& ref, const varCoord& q)
{
	if(ref.pos != q.pos)
		return 0;
	if(ref.a1 == q.a1 && ref.a2 == q.a2)
		return 3;
	else if(ref.a1 == q.a2 && ref.a2 == q.a1)
		return 1;
	return 0;
}

inline int stat_con(const varCoord& ref, const varCoord& q)
{
	if(ref.pos != q.pos)
		return 0;
	if(ref.a1 == q.a1 && ref.a2 == q.a2)
		return 3;
	return 0;
}

static int (*stat_func)(const varCoord& ref, const varCoord& q) = stat_amb;

class snp
{
public:
	snp()
	{
		idxs = new int[npop * 3];
		memset(idxs, 0, sizeof(int) * 3 * npop);
		stats = new double[npop * 4];
		memset(stats, 0, sizeof(double) * 4 * npop);
	}
	~snp()
	{
		delete []idxs;
		delete []stats;
	}

	snp(const snp& dat)
	{
		if(this != & dat)
		{
			coord = dat.coord;
			idxs = new int[npop * 3];
			memcpy(idxs, dat.idxs, sizeof(int) * npop * 3);
			stats = new double[npop * 4];
			memcpy(stats, dat.stats, sizeof(double) * npop * 4);
		}
	}

	snp & operator=(const snp& dat)
	{
		if(this != & dat)
		{
			coord = dat.coord;
			memcpy(idxs, dat.idxs, sizeof(int) * npop * 3);
			memcpy(stats, dat.stats, sizeof(double) * npop * 4);
		}
		return *this;
	}

	void ref_set(int popIdx, const ldref& s)
	{
		coord = s.coord;
		idxs[popIdx * 3] = 3;
		idxs[popIdx * 3 + 1] = s.idx;
		stats[popIdx * 4] = 1 - s.frq;
	}
	void set(int popIdx, const ldref& s)
	{
		idxs[popIdx * 3] = stat_func(coord, s.coord);
		idxs[popIdx * 3 + 1] = s.idx;
		stats[popIdx * 4] = 1 - s.frq;
	}
	void set_beta(int popIdx, const sumstats& s)
	{
		int st = stat_func(coord, s.coord);

		if(st)
		{
			if(st == idxs[popIdx * 3])
				stats[popIdx * 4 + 1] = s.beta;
			else
				stats[popIdx * 4 + 1] = - s.beta;
			idxs[popIdx * 3] += 4;
			idxs[popIdx * 3 + 2] = s.idx;
			stats[popIdx * 4 + 2] = s.pval;
			//stats[popIdx * 4 + 3] = - log10(s.pval);
			stats[popIdx * 4 + 3] = s.log10p;
		}
	}
	bool var() const
	{
		for(int i = 0 ; i < npop; ++i)
			if((idxs[i * 3] & 5) == 5)
				return true;
		return false;
	}
	bool var(int idxPop) const
	{
		return ((idxs[idxPop * 3] & 5) == 5);
	}
	/*
	int ldflag(int idxPop) const
	{
		return 1;
		if((idxs[idxPop * 3] & 3) == 3)
			return 1;
		else if((idxs[idxPop * 3] & 3) == 1)
			return -1;
		else
			return 0;
	}
	*/

	varCoord coord;

	//[[1 in LD ref file, 2 non-rev in LD ref file, 4 in gwas file]
	//[index of LD ref file] [index of gwas file]]
	int *idxs;

	//[[frq] [beta] [p] [logp]]
	double *stats;

	static softpar par;
};

inline std::ostream & operator << (std::ostream & os, const snp& var)
{
	//SNP BP
	os << var.coord.id << '\t' << var.coord.pos;

	//REF_ALLELE
	if((var.idxs[0] & 7) == 7)
		os << '\t' << var.coord.a1;
	else if((var.idxs[0] & 7) == 5)
		os << '\t' << var.coord.a2;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 7) == 7)
			os << ',' << var.coord.a1;
		else if((var.idxs[i * 3] & 7) == 5)
			os << ',' << var.coord.a2;
		else
			os << ",NA";

	//ALT_ALLELE
	if((var.idxs[0] & 7) == 7)
		os << '\t' << var.coord.a2;
	else if((var.idxs[0] & 7) == 5)
		os << '\t' << var.coord.a1;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 7) == 7)
			os << ',' << var.coord.a2;
		else if((var.idxs[i * 3] & 7) == 5)
			os << ',' << var.coord.a1;
		else
			os << ",NA";

	//REF_FRQ
	if((var.idxs[0] & 5) == 5)
		os << '\t' << var.stats[0];
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << var.stats[i * 4];
		else
			os << ",NA";

	//BETA
	double sqrtfrq(sqrt(2 * var.stats[0] * (1 - var.stats[0])));
	if((var.idxs[0] & 5) == 5)
		os << '\t' << var.stats[1] / sqrtfrq;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
	{
		sqrtfrq = sqrt(2 * var.stats[i * 4] * (1 - var.stats[i * 4]));
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << var.stats[i * 4 + 1] / sqrtfrq;
		else
			os << ",NA";
	}

	//SE
	sqrtfrq = sqrt(2 * var.stats[0] * (1 - var.stats[0]));
	if((var.idxs[0] & 5) == 5)
		os << '\t' << 1.0 / sqrt(snp::par.n_gwas[0]) / sqrtfrq;
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
	{
		sqrtfrq = sqrt(2 * var.stats[i * 4] * (1 - var.stats[i * 4]));
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << 1.0 / sqrt(snp::par.n_gwas[i]) / sqrtfrq;
		else
			os << ",NA";
	}

	//-LOG10P
	if((var.idxs[0] & 5) == 5)
		os << '\t' << var.stats[3];
	else
		os << "\tNA";
	for(int i = 1 ; i < npop; ++i)
		if((var.idxs[i * 3] & 5) == 5)
			os << ',' << var.stats[i * 4 + 3];
		else
			os << ",NA";

    return os;
}

class dataset
{
public:
	dataset() : beta(NULL), pval(NULL), mkIdx(NULL), ind(NULL), nsnp(0)
	{
		tau_sq = new double[npop];
		ld = new LDTYPE**[npop]();
	}

	~dataset()
	{
		if(npop)
		{
			if(nsnp)
			{
				for(int i = 0 ; i < npop; ++i)
				{
					if(ld[i] == nullptr) continue;
					for(int j = 0 ; j < nsnp; ++j)
					{
						if(ld[i][j] != nullptr)
							delete []ld[i][j];
					}
				}
			}
			if(tau_sq) delete []tau_sq;
			if(beta) delete []beta;
			if(pval) delete []pval;
			if(mkIdx) delete []mkIdx;
			if(ind) delete []ind;
			if(ld) delete []ld;
		}
	}

	void load(const softpar& par);

	// Load dataset from in-memory arrays. Caller must set the global `npop` to the number of populations
	// before constructing the dataset instance (the constructor uses `npop` for initial allocations).
	void load_from_memory(const softpar& par,
					int _nsnp,
					const double* beta_in,    // length npop * nsnp, row-major (pop major)
					const double* pval_in,    // length npop * nsnp
					const char* ind_in,       // length npop * nsnp
					LDTYPE*** ld_in,         // [npop][nsnp][nsnp] dense matrices
					const int* mkIdx_in);    // length nsnp

	//pop
	double *tau_sq;

	//pop X snp
	double *beta, *pval;
	//snp
	int *mkIdx;
	//pop X snp
	char *ind;

	//pop snp snp
	LDTYPE ***ld;
	std::vector<snp> mks;

	int nsnp;
};


#endif /* DATA_HPP_ */
