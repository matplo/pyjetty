#ifndef PYJETTY_ALEPHR_HH
#define PYJETTY_ALEPHR_HH

#include <string>
#include <fstream>
#include <vector>

#include <TObject.h>

namespace AlephR
{
	int dump(const char *fname, int nevents, bool noparts = false);

	class Particle
	{
	public:
		Particle();
		Particle(const Particle &p);
		Particle(double px, double py, double pz, double m, double q, int pwflag, double d0, double z0, int ntpc, int nitc, int nvdet);

		double 	px() 	 const {return fpx;}
		double 	py() 	 const {return fpy;}
		double 	pz() 	 const {return fpz;}
		double 	m() 	 const {return fm;}
		double 	q() 	 const {return fq;}
		int 	pwflag() const {return fpwflag;}
		double 	d0() 	 const {return fd0;}
		double 	z0() 	 const {return fz0;}
		int 	ntpc() 	 const {return fntpc;}
		int 	nitc() 	 const {return fnitc;}
		int 	nvdet()  const {return fnvdet;}
		double  e()      const {return fe;}
		double  E()      const {return fe;}
		double  pt()     const {return fpt;}

		~Particle() {;}

		void dump() const;
		double calc_e() const;
		double calc_pt() const;

		const std::vector<double> as_vector() const;
		static std::vector<std::string> descr();

	private:
		double 	fpx;
		double 	fpy;
		double 	fpz;
		double 	fm;
		double 	fq;
		int    	fpwflag;
		double 	fd0;
		double 	fz0;
		int 	fntpc;
		int 	fnitc;
		int 	fnvdet;

		double  fe;
		double  fpt;
	};

	class EventHeader
	{
	public:
		EventHeader();
		EventHeader(const EventHeader &h);
		EventHeader(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey);

		void reset(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey);

		int 	run() 	const {return frun;}
		int 	n() 	const {return fn;}
		double 	e() 	const {return fe;}
		int 	vflag() const {return fvflag;}
		double  vx() 	const {return fvx;}
		double  vy() 	const {return fvy;}
		double  ex() 	const {return fex;}
		double  ey() 	const {return fey;}

		void clear();
		void dump() const;

		~EventHeader() {;}

	private:
		int 	frun;
		int 	fn;
		double 	fe;
		int 	fvflag;
		double  fvx;
		double  fvy;
		double  fex;
		double  fey;
	};

	class Event
	{
	public:
		Event();
		Event(const Event &e);
		Event(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey);

		EventHeader get_header() const;

		std::vector<Particle> get_particles() const;
		std::vector< std::vector<double> > get_particles_vdoubles() const;

		void add_particle(const Particle &p);
		void add_particle(double px, double py, double pz, double m, double q, int pwflag, double d0, double z0, int ntpc, int nitc, int nvdet);

		void reset(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey);
		void clear();
		void dump(bool noparts = false) const;

		~Event() {;}

	private:
		EventHeader fHeader;
		std::vector<Particle> fparticles;
	};

	class Reader
	{
	public:
		Reader();
		Reader(const char *fname);
		void open(const char *fname);
		void reset();
		bool read_next_event();
		const Event& get_event();
		~Reader();
	private:
		std::string 	fName;
		std::ifstream 	fStream;
		Event 			fEvent;
	};

	class ReaderLines
	{
	public:
		ReaderLines(const std::vector<std::string> &data);
		bool read_next_event();
		const Event& get_event();
		~ReaderLines();
	private:
		ReaderLines();
		std::vector<std::string> fdata;
		std::vector<std::string>::iterator iter;
		Event 			fEvent;
	};

	void write_root_tree_lines(const std::vector<std::string> &data, const char *outputfname, Size_t nev = -1);
};

#endif
