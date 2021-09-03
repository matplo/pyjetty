// ALEPH_DATA RUN = 14209 EVENT  4218 ECM =   91.276 GEV
//  Primary vertex info flag =  4 vx = -0.0855 vy =  0.0485 ex =  0.0057 ey =  0.0000
// px=    0.001 py=    0.150 pz=   -0.126 m=    0.140 charge -1.0 pwflag    0 d0   0.051 z0  -0.359 ntpc   4 nitc   7 nvdet   2
// px=    3.001 py=   -2.409 pz=   -2.398 m=    0.140 charge  1.0 pwflag    0 d0  -0.013 z0  -0.393 ntpc  21 nitc   5 nvdet   2
// px=    4.377 py=   -3.291 pz=   -3.063 m=    0.140 charge -1.0 pwflag    0 d0   0.000 z0  -0.398 ntpc  18 nitc   3 nvdet   2
// px=    3.602 py=   -2.626 pz=   -2.855 m=    0.140 charge  1.0 pwflag    0 d0   0.001 z0  -0.394 ntpc  20 nitc   3 nvdet   2
// px=   -3.368 py=    2.556 pz=    1.994 m=    0.140 charge -1.0 pwflag    0 d0  -0.047 z0  -0.408 ntpc  21 nitc   0 nvdet   2
// px=    0.978 py=   -0.599 pz=   -0.757 m=    0.140 charge -1.0 pwflag    0 d0   0.000 z0  -0.357 ntpc  20 nitc   5 nvdet   2
// px=   -0.897 py=    0.522 pz=    0.484 m=    0.140 charge -1.0 pwflag    0 d0   0.006 z0  -0.396 ntpc  19 nitc   7 nvdet   2
// px=   -5.595 py=    3.948 pz=    3.960 m=    0.140 charge  1.0 pwflag    0 d0   0.396 z0  -0.340 ntpc  12 nitc   2 nvdet   2
// px=    0.755 py=   -0.068 pz=   -0.025 m=    0.140 charge -1.0 pwflag    0 d0  -0.005 z0  -0.417 ntpc  19 nitc   6 nvdet   2
// px=    5.527 py=   -3.759 pz=   -3.514 m=    0.140 charge -1.0 pwflag    0 d0   0.000 z0  -0.399 ntpc  19 nitc   4 nvdet   2
// px=   -0.324 py=    0.067 pz=    0.457 m=    0.140 charge  1.0 pwflag    0 d0  -0.036 z0  -0.363 ntpc  13 nitc   8 nvdet   0
// px=    7.434 py=   -5.007 pz=   -5.154 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=    0.831 py=   -0.580 pz=   -0.457 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=   -0.658 py=    0.505 pz=    0.500 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=   -2.131 py=    1.651 pz=    1.799 m=    0.025 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=   -2.852 py=    2.115 pz=    2.425 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=   -7.184 py=    4.718 pz=    4.988 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=   -0.602 py=    0.368 pz=    0.525 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=    1.634 py=   -1.048 pz=   -1.371 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=    0.353 py=   -0.387 pz=   -0.731 m=    0.000 charge  0.0 pwflag    4 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=   -8.746 py=    7.889 pz=    7.007 m=    0.052 charge  0.0 pwflag    5 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=    1.140 py=   -0.210 pz=   -0.254 m=    0.015 charge  0.0 pwflag    5 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// px=    0.350 py=   -0.384 pz=   -0.726 m=    0.000 charge  0.0 pwflag    5 d0  -1.000 z0  -1.000 ntpc  -1 nitc  -1 nvdet  -1
// END_EVENT

#include "strutil.hh"
#include "aleph.hh"

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace Aleph
{
	Particle::Particle()
			: fpx(0), fpy(0), fpz(0), fm(0), fq(0), fpwflag(0), fd0(0), fz0(0), fntpc(0), fnitc(0), fnvdet(0), fe(calc_e()), fpt(calc_pt())
	{
		;
	}

	Particle::Particle(const Particle &p)
			: fpx(p.fpx), fpy(p.fpy), fpz(p.fpz), fm(p.fm), fq(p.fq), fpwflag(p.fpwflag), fd0(p.fd0), fz0(p.fz0), fntpc(p.fntpc), fnitc(p.fnitc), fnvdet(p.fnvdet), fe(calc_e()), fpt(calc_pt())
	{
		;
	}


	Particle::Particle(double px, double py, double pz, double m, double q, int pwflag, double d0, double z0, int ntpc, int nitc, int nvdet)
			: fpx(px), fpy(py), fpz(pz), fm(m), fq(q), fpwflag(pwflag), fd0(d0), fz0(z0), fntpc(ntpc), fnitc(nitc), fnvdet(nvdet), fe(calc_e()), fpt(calc_pt())
	{
		;
	}

	void Particle::dump() const
	{
		std::cout << std::showpoint << std::setprecision(2)
			<< std::setw(7) << fpt << " " << std::setw(7) << fe << " " << std::setw(7) << fpx << " " << std::setw(7) << fpy << " " << std::setw(7) << fpz << " " << std::setw(7) << fm << " " << std::setw(7) << fq << " " << std::setw(7) << fpwflag << " " << std::setw(7) << fd0 << " " << std::setw(7) << fz0 << " " << std::setw(7) << fntpc << " " << std::setw(7) << fnitc << " " << std::setw(7) << fnvdet << std::endl;
	}

	double Particle::calc_e() const
	{
		double _fe = std::sqrt(fpx*fpx + fpy*fpy + fpz*fpz + fm*fm);
		return _fe;
	}

	double Particle::calc_pt() const
	{
		double _fpt = std::sqrt(fpx*fpx + fpy*fpy);
		return _fpt;
	}

	const std::vector<double> Particle::as_vector() const
	{
		const double _tmp[] = {	fpx, fpy, fpz, fm, fq, 
								static_cast<double>(fpwflag), 
								fd0, fz0, 
								static_cast<double>(fntpc), 
								static_cast<double>(fnitc), 
								static_cast<double>(fnvdet), 
								fe, fpt};
		std::vector<double> v(_tmp, _tmp + sizeof(_tmp) / sizeof(_tmp[0]));
		return v;
	}

	std::vector<std::string> Particle::descr()
	{
		const std::string _tmp[] = {std::string("px"), std::string("py"), std::string("pz"), std::string("m"), 
									std::string("q"), std::string("pwflag"), std::string("d0"), std::string("z0"), 
									std::string("ntpc"), std::string("nitc"), std::string("nvdet"), 
									std::string("e"), std::string("pt")};
		std::vector<std::string> v(_tmp, _tmp + sizeof(_tmp) / sizeof(_tmp[0]));
		return v;		
	}

	// ----------
	EventHeader::EventHeader()
		: frun(0), fn(0), fe(0), fvflag(0), fvx(0), fvy(0), fex(0), fey(0)
	{
		;
	}

	EventHeader::EventHeader(const EventHeader& h)
		: frun(h.frun), fn(h.fn), fe(h.fe), fvflag(h.fvflag), fvx(h.fvx), fvy(h.fvy), fex(h.fex), fey(h.fey)
	{
		;
	}

	EventHeader::EventHeader(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey)
		: frun(run), fn(n), fe(e), fvflag(vflag), fvx(vx), fvy(vy), fex(ex), fey(ey)
	{
		;
	}

	void EventHeader::dump() const
	{
		std::cout << "EventHeader: " << std::setw(7) << frun << " " << std::setw(7) << fn << " " << std::setw(7) << fe << " " << std::setw(7) << fvflag << " " << std::setw(7) << fvx << " " << std::setw(7) << fvy << " " << std::setw(7) << fex << " " << std::setw(7) << fey << std::endl;
	}


	void EventHeader::clear()
	{
		reset(0, 0, 0, 0, 0, 0, 0, 0);
	}

	void EventHeader::reset(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey)
	{
		frun   = run;
		fn     = n;
		fe     = e;
		fvflag = vflag;
		fvx    = vx;
		fvy    = vy;
		fex    = ex;
		fey    = ey;
	}
	// ----------

	Event::Event()
		: fHeader(), fparticles()
	{
		;
	}

	Event::Event(const Event &e)
		: fHeader(e.fHeader), fparticles(e.fparticles)
	{
		;
	}

	Event::Event(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey)
		: fHeader(run, n, e, vflag, vx, vy, ex, ey), fparticles()
	{
		;
	}

	void Event::dump(bool noparts) const
	{
		fHeader.dump();
		if (noparts == false)
		{
			std::cout << std::setw(7) << "    pt" << std::setw(7) << "    e" << std::setw(7) <<  "  px" << std::setw(7) << "  py" << std::setw(7) << "  pz" << std::setw(7) << "  mass " << std::setw(7) << "  charge" << std::setw(7) << "  pwflag" << std::setw(7) << "    d0" << std::setw(7) << "    z0" << std::setw(7) << "  ntpc" << std::setw(7) << "  nitc" << std::setw(7) << "  nvdet" << std::endl;
			for (unsigned int i = 0; i < fparticles.size(); i++)
				fparticles[i].dump();
			std::cout << "---" << std::endl;
		}
	}

	std::vector<Particle> Event::get_particles() const {return fparticles;}
	std::vector< std::vector<double> > Event::get_particles_vdoubles() const 
	{	
		std::vector< std::vector<double> > vv;
		for (const Particle &p : fparticles)
		{
			const std::vector<double> &v = p.as_vector();
			vv.push_back(v);
		}
		return vv;
	}

	EventHeader Event::get_header() const {return fHeader;}

	void Event::add_particle(const Particle &p)
	{
		fparticles.push_back(p);
	}

	void Event::add_particle(double px, double py, double pz, double m, double q, int pwflag, double d0, double z0, int ntpc, int nitc, int nvdet)
	{
		Particle p(px, py, pz, m, q, pwflag, d0, z0, ntpc, nitc, nvdet);
		fparticles.push_back(p);
	}

	void Event::clear()
	{
		fHeader.clear();
		fparticles.clear();
	}

	void Event::reset(int run, int n, double e, int vflag, double vx, double vy, double ex, double ey)
	{
		fHeader.reset(run, n, e, vflag, vx, vy, ex, ey);
		fparticles.clear();
	}

	// ----------

	Reader::Reader()
		: fName(), fStream(), fEvent()
	{
		;
	}

	Reader::~Reader()
	{
		fEvent.clear();
		fStream.close();
	}

	Reader::Reader(const char *fname)
		: fName(fname), fStream(), fEvent()
	{
		fStream.open(fName);
	}

	const Event& Reader::get_event() {return fEvent;}

	void Reader::reset()
	{
		fEvent.clear();
		fStream.close();
	}

	void Reader::open(const char *fname)
	{
		fStream.close();
		fName = fname;
		fStream.open(fName);
	}

	bool Reader::read_next_event()
	{
		if (!fStream.good())
			return false;

		std::string line;

		bool reading_parts = false;
		bool reading_event = false;

		int run    = -99;
		int event  = -99;
		double ecm = -99.0;
		int vflag  = -99;
		double vx  = -99.0;
		double vy  = -99.0;
		double ex  = -99.0;
		double ey  = -99.0;

		while (std::getline(fStream, line))
		{
			if (line.find("ALEPH_DATA", 0) != std::string::npos)
			{
				// std::cout << line << " - start of the event found" << std::endl;
				if (reading_event)
				{
					std::cerr << "Error: something went wrong : new event before finishing the last one" << event << std::endl;
				}
				reading_event = true;
				StrUtil::replace_substring(line, " ", "");
				StrUtil::replace_substring(line, "EVENT", " EVENT=");
				StrUtil::replace_substring(line, "ECM", " ECM");
				StrUtil::replace_substring(line, "GEV", "");

				StrUtil::Args _args(line);
				run   = _args.getI("ALEPH_DATARUN", 0);
				event = _args.getI("EVENT", 0);
				ecm   = _args.getD("ECM", 0.0);

				fEvent.clear();
			}
			else
			{
				if (reading_event)
				{
					if (line.find("Primary vertex info flag", 0) != std::string::npos)
					{
						StrUtil::replace_substring(line, "Primary vertex info flag", "vflag");
						StrUtil::replace_substring(line, " ", "");
						StrUtil::replace_substring(line, "vx", " vx");
						StrUtil::replace_substring(line, "vy", " vy");
						StrUtil::replace_substring(line, "ex", " ex");
						StrUtil::replace_substring(line, "ey", " ey");
						reading_parts = true;
						// std::cout << line << std::endl;
						StrUtil::Args _args(line);
						vflag = _args.getI("vflag", -99);
						vx    = _args.getD("vx", -99.);
						vy    = _args.getD("vy", -99.);
						ex    = _args.getD("ex", -99.);
						ey    = _args.getD("ey", -99.);
						// std::cout << vflag << " " << std::setw(7) << vx << " " << std::setw(7) << vy << " " << std::setw(7) << ex << " " << std::setw(7) << ey << std::endl;
					}
					else if (reading_parts)
					{
						if (line.find("px=", 0) == 0)
						{
							if (fEvent.get_particles().size() == 0)
							{
								fEvent.reset(run, event, ecm, vflag, vx, vy, ex, ey);
							}
							StrUtil::replace_substring(line, " ", "");
							std::vector<std::string> vars = StrUtil::split_to_vector("px=,py=,pz=,m=,charge,pwflag,d0,z0,ntpc,nitc,nvdet", ",");
							std::vector<std::string> tovars = StrUtil::split_to_vector("px=, py=, pz=, m=, charge=, pwflag=, d0=, z0=, ntpc=, nitc=, nvdet=", ",");
							for ( unsigned int i = 0; i < vars.size(); i++)
							{
								StrUtil::replace_substring(line, vars[i].c_str(), tovars[i].c_str());
							}
							StrUtil::Args _args(line);
							double px     = _args.getD("px", 0);
							double py     = _args.getD("py", 0);
							double pz     = _args.getD("pz", 0);
							double m      = _args.getD("m", 0);
							double q      = _args.getD("charge", 0);
							   int pwflag = _args.getI("pwflag", 0);
							double d0     = _args.getD("d0", 0);
							double z0     = _args.getD("z0", 0);
							   int ntpc   = _args.getI("ntpc", 0);
							   int nitc   = _args.getI("nitc", 0);
							   int nvdet  = _args.getI("nvdet", 0);
							Particle p(px, py, pz, m, q, pwflag, d0, z0, ntpc, nitc, nvdet);
							fEvent.add_particle(p);
						}
						else
						{
							if (line.find("END_EVENT", 0) != std::string::npos)
							{
								// this is ok
								// std::cout << line << " - end of the event found - while reading particles" << std::endl;
							}
							else
							{
								std::cerr << "Error: while reading particles - unknown tag found" << std::endl;
								std::cerr << "Error: - line is: " << line << std::endl;
							}
						}
					}
				}
			}
			// std::cout << line << std::endl;
			if (line.find("END_EVENT", 0) != std::string::npos)
			{
				// std::cout << line << " - end of the event found" << std::endl;
				if (reading_event && reading_parts)
				{
					reading_event = false;
					reading_parts = false;
					return true;
				}
				else
				{
					std::cerr << "Error: Something awfully wrong with the event structure..." << std::endl;
					return false;
					break;
				}
			}
		}
		fStream.close();
		return false;
	}

	int dump(const char *fname, int nevents, bool noparts)
	{
		int n = 0;
		std::cout << "Dumping file: " << fname << std::endl;
		Reader r(fname);
		while (r.read_next_event())
		{
			const Event &ev = r.get_event();
			ev.dump(noparts);
			n++;
			if (nevents > 0)
				if (n >= nevents)
					break;
		}
		return 0;
	}

	// ----------

	ReaderLines::ReaderLines()
		: fdata()
		, iter(fdata.begin())
	{
		;
	}

	ReaderLines::~ReaderLines()
	{
		;
	}

	ReaderLines::ReaderLines(const std::vector<std::string> &data)
		: fdata(data)
		, iter(fdata.begin())
	{
		;
	}

	const Event& ReaderLines::get_event() {return fEvent;}

	bool ReaderLines::read_next_event()
	{
		std::string line;

		bool reading_parts = false;
		bool reading_event = false;

		int run    = -99;
		int event  = -99;
		double ecm = -99.0;
		int vflag  = -99;
		double vx  = -99.0;
		double vy  = -99.0;
		double ex  = -99.0;
		double ey  = -99.0;

		// while (auto itline = next(iter))
		for (std::vector<std::string>::iterator it = iter; it != fdata.end(); ++it)
		{
			iter = it;
			std::string line = *it;
			if (line.find("ALEPH_DATA", 0) != std::string::npos)
			{
				// std::cout << line << " - start of the event found" << std::endl;
				if (reading_event)
				{
					std::cerr << "Error: something went wrong : new event before finishing the last one" << event << std::endl;
				}
				reading_event = true;
				StrUtil::replace_substring(line, " ", "");
				StrUtil::replace_substring(line, "EVENT", " EVENT=");
				StrUtil::replace_substring(line, "ECM", " ECM");
				StrUtil::replace_substring(line, "GEV", "");

				StrUtil::Args _args(line);
				run   = _args.getI("ALEPH_DATARUN", 0);
				event = _args.getI("EVENT", 0);
				ecm   = _args.getD("ECM", 0.0);

				fEvent.clear();
			}
			else
			{
				if (reading_event)
				{
					if (line.find("Primary vertex info flag", 0) != std::string::npos)
					{
						StrUtil::replace_substring(line, "Primary vertex info flag", "vflag");
						StrUtil::replace_substring(line, " ", "");
						StrUtil::replace_substring(line, "vx", " vx");
						StrUtil::replace_substring(line, "vy", " vy");
						StrUtil::replace_substring(line, "ex", " ex");
						StrUtil::replace_substring(line, "ey", " ey");
						reading_parts = true;
						// std::cout << line << std::endl;
						StrUtil::Args _args(line);
						vflag = _args.getI("vflag", -99);
						vx    = _args.getD("vx", -99.);
						vy    = _args.getD("vy", -99.);
						ex    = _args.getD("ex", -99.);
						ey    = _args.getD("ey", -99.);
						// std::cout << vflag << " " << std::setw(7) << vx << " " << std::setw(7) << vy << " " << std::setw(7) << ex << " " << std::setw(7) << ey << std::endl;
					}
					else if (reading_parts)
					{
						if (line.find("px=", 0) == 0)
						{
							if (fEvent.get_particles().size() == 0)
							{
								fEvent.reset(run, event, ecm, vflag, vx, vy, ex, ey);
							}
							StrUtil::replace_substring(line, " ", "");
							std::vector<std::string> vars = StrUtil::split_to_vector("px=,py=,pz=,m=,charge,pwflag,d0,z0,ntpc,nitc,nvdet", ",");
							std::vector<std::string> tovars = StrUtil::split_to_vector("px=, py=, pz=, m=, charge=, pwflag=, d0=, z0=, ntpc=, nitc=, nvdet=", ",");
							for ( unsigned int i = 0; i < vars.size(); i++)
							{
								StrUtil::replace_substring(line, vars[i].c_str(), tovars[i].c_str());
							}
							StrUtil::Args _args(line);
							double px     = _args.getD("px", 0);
							double py     = _args.getD("py", 0);
							double pz     = _args.getD("pz", 0);
							double m      = _args.getD("m", 0);
							double q      = _args.getD("charge", 0);
							   int pwflag = _args.getI("pwflag", 0);
							double d0     = _args.getD("d0", 0);
							double z0     = _args.getD("z0", 0);
							   int ntpc   = _args.getI("ntpc", 0);
							   int nitc   = _args.getI("nitc", 0);
							   int nvdet  = _args.getI("nvdet", 0);
							Particle p(px, py, pz, m, q, pwflag, d0, z0, ntpc, nitc, nvdet);
							fEvent.add_particle(p);
						}
						else
						{
							if (line.find("END_EVENT", 0) != std::string::npos)
							{
								// this is ok
								// std::cout << line << " - end of the event found - while reading particles" << std::endl;
								iter++;
							}
							else
							{
								std::cerr << "Error: while reading particles - unknown tag found" << std::endl;
								std::cerr << "Error: - line is: " << line << std::endl;
							}
						}
					}
				}
			}
			// std::cout << line << std::endl;
			if (line.find("END_EVENT", 0) != std::string::npos)
			{
				// std::cout << line << " - end of the event found" << std::endl;
				if (reading_event && reading_parts)
				{
					reading_event = false;
					reading_parts = false;
					return true;
				}
				else
				{
					std::cerr << "Error: Something awfully wrong with the event structure..." << std::endl;
					return false;
					break;
				}
			}
		}
		return false;
	}
};
